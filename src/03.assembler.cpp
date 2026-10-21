#include <hip/debruijn_core.hpp>
#include <hip/panel_core.hpp>
#include <hip/validator_core.hpp>

#include <parasail.h>

#include <tsl/robin_set.h>
#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

namespace fs = std::filesystem;

struct Args {
    std::string input_dir;
    std::string mers_file;
    std::string output = "assembled.fasta";
    int kmer = 19;
    // Reads-mode default: single k=27. `--multi-k` opts into multi-k assembly;
    // `-k` or `--no-multi-k` forces single-k mode.
    std::vector<int> multi_k = {27};
    int min_count = 1;
    int max_reads = 200;
    int max_loci_gap = 5000;
    int max_assign = 3;   // load-time cap; 0 = keep all
    int threads = 4;
    tsl::robin_set<int> debug_eles;
    bool trim = true;
    int trim_margin = 10;   // bp kept beyond the first/last element-specific k-mer
    bool polish = true;   // reads-mode post-assembly consensus polish
    bool trim_set = false;
    bool no_trim_set = false;
    double repeat_ratio = 1.2;
    bool output_snp = false;
    bool help_requested = false;
    std::string cne_ref_path;
    // Rescue/panel path (reads mode): score candidate contigs and low-support
    // reads against the element's leave-one-out panels. Disabled when empty.
    std::string panels_path;
    double rescue_tau = 0.65;
    int rescue_min_distinct = 2;
    int rescue_min_span = 30;
    int min_output_distinct = 12;
    int min_fallback_distinct = 16;
    bool rescue_lowconf = true;
    int rescue_min_support = 3;
    double min_rescue_score = 1.0;    // 0 = keep all T3 rescue-tagged reads
    bool rescue_enabled = true;
};

struct Read {
    std::string seq_id;
    std::string seq;
    int strand;
    int ele_id;
    std::string chr;
    int start = 0;
    int end = 0;
    double rescue = 0.0;   // T3 rescue score from validate (0 = strict output)
};

// Uppercase, N-mask and (for the minus strand) reverse-complement a validated
// sequence so every downstream consumer sees element orientation.
static inline void normalize_read_seq(std::string& seq, int strand) {
    for (char& c : seq) {
        if (c >= 'a' && c <= 'z') c = static_cast<char>(c - 'a' + 'A');
    }
    for (char& c : seq) {
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T') c = 'N';
    }
    if (strand == -1) {
        thread_local std::string rc;
        reverse_complement(seq, rc);
        seq = rc;
    }
}

// Shared output evidence floors (both assembly modes):
//   - a de Bruijn contig needs >= 12 distinct element mers, or >= 8 when the
//     contig verifies against the cross-species CNE reference;
//   - a single-sequence claim (reads-mode fallback, genome single-window
//     cluster) needs >= 16 distinct element mers on its own.
static constexpr int kMinOutputDistinct = 12;
static constexpr int kMinOutputDistinctLow = 8;
static constexpr size_t kMinOutputSupportLow = 4;
static constexpr double kMinCneRefIdentity = 0.75;
static constexpr int kMinFallbackDistinct = 16;

static inline bool asm_debug_enabled() {
    static const bool enabled = std::getenv("MERCURY_ASM_DEBUG") != nullptr;
    return enabled;
}

static inline bool rescue_handoff_enabled() {
    static const bool enabled = std::getenv("MERCURY_NO_RESCUE_HANDOFF") == nullptr;
    return enabled;
}

static inline double min_cne_ref_identity() {
    static const double v = []() {
        const char* s = std::getenv("MERCURY_CNE_IDENT");
        return s ? std::atof(s) : kMinCneRefIdentity;
    }();
    return v;
}

Args parse_args(int argc, char* argv[]) {
    Args args;
    int i = 1;
    while (i < argc) {
        std::string arg(argv[i]);
        if (arg == "-h" || arg == "--help") {
            args.help_requested = true;
            return args;
        } else if (arg == "--mers") {
            if (++i >= argc) throw std::runtime_error("Missing value for --mers");
            args.mers_file = argv[i];
        } else if (arg == "-o" || arg == "--output") {
            if (++i >= argc) throw std::runtime_error("Missing value for --output");
            args.output = argv[i];
        } else if (arg == "-k" || arg == "--kmer") {
            if (++i >= argc) throw std::runtime_error("Missing value for --kmer");
            args.kmer = std::stoi(argv[i]);
            args.multi_k.clear();   // explicit -k means single-k mode
        } else if (arg == "--no-multi-k") {
            args.multi_k.clear();
        } else if (arg == "--multi-k") {
            if (++i >= argc) throw std::runtime_error("Missing value for --multi-k");
            std::string list = argv[i];
            size_t p = 0;
            while (p < list.size()) {
                size_t q = list.find(',', p);
                if (q == std::string::npos) q = list.size();
                if (q > p) {
                    int k = std::stoi(list.substr(p, q - p));
                    if (k < 3) throw std::runtime_error("--multi-k values must be >= 3");
                    args.multi_k.push_back(k);
                }
                p = q + 1;
            }
        } else if (arg == "--min-count") {
            if (++i >= argc) throw std::runtime_error("Missing value for --min-count");
            args.min_count = std::stoi(argv[i]);
        } else if (arg == "--max-reads") {
            if (++i >= argc) throw std::runtime_error("Missing value for --max-reads");
            args.max_reads = std::stoi(argv[i]);
        } else if (arg == "--max-loci-gap") {
            if (++i >= argc) throw std::runtime_error("Missing value for --max-loci-gap");
            args.max_loci_gap = std::stoi(argv[i]);
        } else if (arg == "--max-assign") {
            if (++i >= argc) throw std::runtime_error("Missing value for --max-assign");
            args.max_assign = std::stoi(argv[i]);
        } else if (arg == "-t" || arg == "--threads") {
            if (++i >= argc) throw std::runtime_error("Missing value for --threads");
            args.threads = std::stoi(argv[i]);
        } else if (arg == "--debug-ele") {
            if (++i >= argc) throw std::runtime_error("Missing value for --debug-ele");
            std::string list = argv[i];
            size_t p = 0;
            while (p < list.size()) {
                size_t q = list.find(',', p);
                if (q == std::string::npos) q = list.size();
                if (q > p) args.debug_eles.insert(std::stoi(list.substr(p, q - p)));
                p = q + 1;
            }
        } else if (arg == "--trim") {
            args.trim = true;
            args.trim_set = true;
        } else if (arg == "--no-trim") {
            args.trim = false;
            args.no_trim_set = true;
        } else if (arg == "--trim-margin") {
            if (++i >= argc) throw std::runtime_error("Missing value for --trim-margin");
            args.trim_margin = std::stoi(argv[i]);
        } else if (arg == "--no-polish") {
            args.polish = false;
        } else if (arg == "--repeat-ratio") {
            if (++i >= argc) throw std::runtime_error("Missing value for --repeat-ratio");
            args.repeat_ratio = std::stod(argv[i]);
        } else if (arg == "--snp") {
            args.output_snp = true;
        } else if (arg == "--cne-ref") {
            if (++i >= argc) throw std::runtime_error("Missing value for --cne-ref");
            args.cne_ref_path = argv[i];
        } else if (arg == "--panels") {
            if (++i >= argc) throw std::runtime_error("Missing value for --panels");
            args.panels_path = argv[i];
        } else if (arg == "--rescue-tau") {
            if (++i >= argc) throw std::runtime_error("Missing value for --rescue-tau");
            args.rescue_tau = std::stod(argv[i]);
        } else if (arg == "--rescue-min-distinct") {
            if (++i >= argc) throw std::runtime_error("Missing value for --rescue-min-distinct");
            args.rescue_min_distinct = std::stoi(argv[i]);
        } else if (arg == "--rescue-min-span") {
            if (++i >= argc) throw std::runtime_error("Missing value for --rescue-min-span");
            args.rescue_min_span = std::stoi(argv[i]);
        } else if (arg == "--min-output-distinct") {
            if (++i >= argc) throw std::runtime_error("Missing value for --min-output-distinct");
            args.min_output_distinct = std::stoi(argv[i]);
        } else if (arg == "--min-fallback-distinct") {
            if (++i >= argc) throw std::runtime_error("Missing value for --min-fallback-distinct");
            args.min_fallback_distinct = std::stoi(argv[i]);
        } else if (arg == "--rescue-lowconf") {
            args.rescue_lowconf = true;
        } else if (arg == "--rescue-min-support") {
            if (++i >= argc) throw std::runtime_error("Missing value for --rescue-min-support");
            args.rescue_min_support = std::stoi(argv[i]);
        } else if (arg == "--min-rescue-score") {
            if (++i >= argc) throw std::runtime_error("Missing value for --min-rescue-score");
            args.min_rescue_score = std::stod(argv[i]);
        } else if (arg == "--no-rescue") {
            args.rescue_enabled = false;
        } else if (arg[0] != '-') {
            if (args.input_dir.empty()) args.input_dir = arg;
            else throw std::runtime_error("Unexpected argument: " + arg);
        } else {
            throw std::runtime_error("Unknown option: " + arg);
        }
        ++i;
    }
    if (args.trim_set && args.no_trim_set) {
        throw std::runtime_error("Conflicting flags: --trim and --no-trim");
    }
    return args;
}

static inline bool is_gz_file(const std::string& path) {
    return path.size() >= 3 && path.compare(path.size() - 3, 3, ".gz") == 0;
}

class BufferedLineReader {
public:
    explicit BufferedLineReader(const std::string& path) {
        if (is_gz_file(path)) {
            std::string cmd = "pigz -dc \"" + path + "\" 2>/dev/null";
            fp_ = popen(cmd.c_str(), "r");
            is_pipe_ = true;
        } else {
            fp_ = fopen(path.c_str(), "r");
            is_pipe_ = false;
        }
        if (!fp_) {
            throw std::runtime_error("Cannot open file: " + path);
        }
    }

    ~BufferedLineReader() {
        if (fp_) {
            if (is_pipe_) pclose(fp_);
            else fclose(fp_);
        }
    }

    bool readline(std::string& line) {
        line.clear();
        while (true) {
            auto nl = buffer_.find('\n', pos_);
            if (nl != std::string::npos) {
                line = buffer_.substr(pos_, nl - pos_);
                pos_ = nl + 1;
                return true;
            }
            if (pos_ < buffer_.size()) {
                buffer_ = buffer_.substr(pos_);
            } else {
                buffer_.clear();
            }
            pos_ = 0;
            size_t bytes_read = fread(raw_, 1, sizeof(raw_), fp_);
            if (bytes_read == 0) {
                if (!buffer_.empty()) {
                    line = std::move(buffer_);
                    buffer_.clear();
                    return true;
                }
                return false;
            }
            buffer_.append(raw_, bytes_read);
        }
    }

private:
    FILE* fp_ = nullptr;
    bool is_pipe_ = false;
    char raw_[65536];
    std::string buffer_;
    size_t pos_ = 0;
};

struct GenomePos {
    std::string chr;
    int start;
    int end;
};

bool parse_genome_id(const std::string& seq_id, GenomePos& pos) {
    const size_t dash = seq_id.rfind('-');
    if (dash == std::string::npos || dash == 0 || dash + 1 >= seq_id.size())
        return false;
    const size_t colon = seq_id.rfind(':', dash);
    if (colon == std::string::npos || colon == 0 || colon + 1 >= dash)
        return false;
    for (size_t i = 0; i < seq_id.size(); ++i) {
        const unsigned char c = static_cast<unsigned char>(seq_id[i]);
        if (c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\f' || c == '\v')
            return false;
    }
    for (size_t i = colon + 1; i < dash; ++i)
        if (seq_id[i] < '0' || seq_id[i] > '9') return false;
    for (size_t i = dash + 1; i < seq_id.size(); ++i)
        if (seq_id[i] < '0' || seq_id[i] > '9') return false;
    pos.chr = seq_id.substr(0, colon);
    pos.start = std::stoi(seq_id.substr(colon + 1, dash - colon - 1));
    pos.end = std::stoi(seq_id.substr(dash + 1));
    return true;
}

struct PanelScore {
    double max_id = -1.0;
    double mean_id = -1.0;
    double median_id = -1.0;
    int n_ge80 = 0;
    int n_ge90 = 0;
    int n_used = 0;
};

struct PanelScorer {
    const std::unordered_map<int, std::vector<std::string>>* panels = nullptr;
    parasail_matrix_t* mat = nullptr;

    bool enabled() const { return panels != nullptr && mat != nullptr && !panels->empty(); }

    PanelScore score(int ele, const std::string& seq, double min_cov = 0.5) const {
        PanelScore out;
        if (!enabled() || seq.empty()) return out;
        auto it = panels->find(ele);
        if (it == panels->end()) return out;
        std::vector<double> ids;
        ids.reserve(it->second.size());
        for (const std::string& ref : it->second) {
            if (ref.empty()) continue;
            parasail_result_t* r = parasail_sw_stats_striped_16(
                seq.c_str(), static_cast<int>(seq.size()),
                ref.c_str(), static_cast<int>(ref.size()), 2, 1, mat);
            if (r->stats->length > 0) {
                const int mn = static_cast<int>(std::min(seq.size(), ref.size()));
                const double cov = mn > 0 ? static_cast<double>(r->stats->length) / mn : 0.0;
                if (cov >= min_cov) {
                    ids.push_back(static_cast<double>(r->stats->matches) / r->stats->length);
                }
            }
            parasail_result_free(r);
        }
        if (ids.empty()) return out;
        std::sort(ids.begin(), ids.end());
        out.n_used = static_cast<int>(ids.size());
        out.max_id = ids.back();
        double sum = 0.0;
        for (double v : ids) {
            sum += v;
            if (v >= 0.80) ++out.n_ge80;
            if (v >= 0.90) ++out.n_ge90;
        }
        out.mean_id = sum / ids.size();
        out.median_id = ids[ids.size() / 2];
        return out;
    }

    double identity(int ele, const std::string& seq, double min_cov = 0.5) const {
        return score(ele, seq, min_cov).max_id;
    }
};

struct ContigScore {
    int hits = 0;
    int distinct = 0;
};

ContigScore score_contig(const std::string& seq, const MerQueryManager& mqm, int ele_id) {
    ContigScore sc;
    int mer_size = mqm.get_mer_size();
    tsl::robin_set<uint32_t> seen;
    for (size_t i = 0; i + mer_size <= seq.size(); ++i) {
        uint32_t code;
        if (!encode_mer_at(seq, i, mer_size, code)) continue;
        auto it = mqm.compressed_mer_query.find(code);
        if (it != mqm.compressed_mer_query.end() && mqm.has_id(it->second, ele_id)) {
            ++sc.hits;
            seen.insert(code);
        }
    }
    sc.distinct = static_cast<int>(seen.size());
    return sc;
}

std::unordered_map<std::string, int> build_node_scores(
    const Graph& graph, const MerQueryManager& mqm, int ele_id, int k)
{
    std::unordered_map<std::string, int> scores;
    scores.reserve(graph.size());
    int mer_size = mqm.get_mer_size();
    for (const auto& [node, _] : graph) {
        if (node.size() < (size_t)mer_size) continue;
        int score = 0;
        for (size_t i = 0; i + mer_size <= node.size(); ++i) {
            uint32_t code;
            if (!encode_mer_at(node, i, mer_size, code)) continue;
            auto it = mqm.compressed_mer_query.find(code);
            if (it != mqm.compressed_mer_query.end() && mqm.has_id(it->second, ele_id)) {
                ++score;
            }
        }
        if (score > 0) {
            scores[node] = score;
        }
    }
    return scores;
}

std::string debruijn_assemble(const std::vector<std::string>& reads, int k,
                               int min_count, const MerQueryManager& mqm, int ele_id,
                               std::ofstream* snp_out = nullptr,
                               std::ofstream* gfa_out = nullptr,
                               bool debug = false)
{
    const bool need_kcounts = (gfa_out != nullptr);
    auto [graph, kmer_counts, e_counts] =
        de_bruijn_graph(reads, k, min_count, need_kcounts, snp_out != nullptr);
    if (graph.empty()) return "";

    std::unordered_map<std::string, int> node_scores = build_node_scores(graph, mqm, ele_id, k);

    std::vector<std::string> contig_path;
    std::string contig = assemble_sequence(graph, node_scores, &contig_path);

    if (!node_scores.empty() && !contig.empty()) {
        ContigScore cs1 = score_contig(contig, mqm, ele_id);
        if (cs1.hits < 3) {
            auto path2 = anchored_trace(graph, node_scores);
            std::string contig2 = path_to_string(path2);
            ContigScore cs2 = score_contig(contig2, mqm, ele_id);
            if (debug)
                std::cerr << "DBG ele=" << ele_id << " pass1 hits=" << cs1.hits
                          << " len=" << contig.size() << " | pass2 hits=" << cs2.hits
                          << " len=" << contig2.size()
                          << (cs2.hits > cs1.hits ? " USE2" : " keep1") << "\n";
            if (cs2.hits > cs1.hits) {
                contig = std::move(contig2);
                contig_path = std::move(path2);
            }
        }
    }

    bool has_var = false;
    if (snp_out && !contig_path.empty())
        has_var = scan_snps(graph, e_counts, contig_path, contig, ele_id, k, *snp_out);

    if (has_var && gfa_out && !contig_path.empty())
        export_path_gfa(graph, node_scores, contig_path, ele_id, k, *gfa_out, &kmer_counts);

    return contig;
}

static int target_span(const std::string& seq, const MerQueryManager& mqm, int ele_id) {
    const int mer_size = mqm.get_mer_size();
    int first = -1, last = -1;
    for (size_t i = 0; i + mer_size <= seq.size(); ++i) {
        uint32_t code;
        if (!encode_mer_at(seq, i, mer_size, code)) continue;
        auto it = mqm.compressed_mer_query.find(code);
        if (it != mqm.compressed_mer_query.end() && mqm.has_id(it->second, ele_id)) {
            if (first < 0) first = static_cast<int>(i);
            last = static_cast<int>(i);
        }
    }
    return first < 0 ? 0 : last - first + mer_size;
}

void trim_contig(std::string& seq, const MerQueryManager& mqm, int ele_id, int margin) {
    if (seq.size() < 26) return;
    int mer_size = mqm.get_mer_size();
    int n = static_cast<int>(seq.size()) - mer_size + 1;

    int first_hit = -1, last_hit = -1;
    for (int i = 0; i < n; ++i) {
        uint32_t code;
        if (!encode_mer_at(seq, static_cast<size_t>(i), mer_size, code)) continue;
        auto it = mqm.compressed_mer_query.find(code);
        if (it != mqm.compressed_mer_query.end() && mqm.has_id(it->second, ele_id)) {
            if (first_hit < 0) first_hit = i;
            last_hit = i;
        }
    }

    if (first_hit < 0) return;

    int start = std::max(0, first_hit - margin);
    int end = std::min(n, last_hit + margin);
    auto trimmed = seq.substr(start, end + mer_size - 1 - start);
    if (trimmed.size() >= 20) {
        seq = std::move(trimmed);
    }
}

// Post-assembly consensus polish. The de Bruijn contig is a single graph path
// (no base voting), so a low-count edge can survive into it. Re-align the
// element's reads to the assembled contig (k-mer-anchored placement, both
// strands) and take a per-position majority vote, iterated. Reads vote only if
// they place at >= min_ident over their overlap, so divergent/paralog reads do
// not corrupt the contig.
static void polish_contig(std::string& contig, const std::vector<std::string>& reads,
                          int iters = 1, double min_ident = 0.95, int k = 13) {
    const int L = static_cast<int>(contig.size());
    if (L < k + 5 || reads.size() < 2) return;
    static const char BASES[4] = {'A', 'C', 'G', 'T'};
    auto bidx = [](char c) -> int {
        switch (c) { case 'A': return 0; case 'C': return 1; case 'G': return 2; case 'T': return 3; default: return -1; }
    };
    auto comp = [](char c) -> char {
        switch (c) { case 'A': return 'T'; case 'C': return 'G'; case 'G': return 'C'; case 'T': return 'A'; default: return 'N'; }
    };
    std::unordered_map<uint32_t, std::vector<int>> idx;
    std::string rcbuf;
    // Flat, epoch-stamped offset counters: offset d lies in [-base, L-k], so a
    // dense array replaces the per-read hash map (no allocation per read).
    thread_local std::vector<int> offcnt, offstamp;
    thread_local std::vector<int> touched;
    thread_local int offep = 0;
    for (int it = 0; it < iters; ++it) {
        idx.clear();
        for (int i = 0; i + k <= L; ++i) {
            uint32_t code;
            if (encode_mer_at(contig, static_cast<size_t>(i), k, code)) idx[code].push_back(i);
        }
        std::vector<std::array<int, 4>> cnt(L);
        for (int i = 0; i < L; ++i) {
            cnt[i] = {0, 0, 0, 0};
            int b = bidx(contig[i]);
            if (b >= 0) cnt[i][b] = 1;   // contig base seeds the vote (tie-break)
        }
        auto place = [&](const std::string& q, int& m, int& ov, int& d) {
            const int n = static_cast<int>(q.size());
            m = -1; ov = 0; d = 0;
            if (n < k) return;
            const int base = n - k;            // offset d_min = -base
            const int D = (L - k) + base + 1;  // offsets in [-base, L-k]
            if (static_cast<int>(offcnt.size()) < D) {
                offcnt.resize(D);
                offstamp.resize(D, 0);
            }
            ++offep;
            if (offep == 0) { std::fill(offstamp.begin(), offstamp.end(), 0); offep = 1; }
            touched.clear();
            for (int i = 0; i + k <= n; ++i) {
                uint32_t code;
                if (!encode_mer_at(q, static_cast<size_t>(i), k, code)) continue;
                auto ii = idx.find(code);
                if (ii == idx.end()) continue;
                for (int cp : ii->second) {
                    int x = cp - i + base;
                    if (offstamp[x] != offep) { offstamp[x] = offep; offcnt[x] = 0; touched.push_back(x); }
                    ++offcnt[x];
                }
            }
            if (touched.empty()) return;
            int bx = touched[0], bc = offcnt[touched[0]];
            for (int x : touched) if (offcnt[x] > bc) { bc = offcnt[x]; bx = x; }
            int bd = bx - base;
            int lo = std::max(0, bd), hi = std::min(L, bd + n);
            ov = hi - lo;
            if (ov <= 0) return;
            int mm = 0;
            for (int p = lo; p < hi; ++p) if (q[p - bd] == contig[p]) ++mm;
            m = mm; d = bd;
        };
        for (const auto& r : reads) {
            rcbuf.assign(r.rbegin(), r.rend());
            for (char& c : rcbuf) c = comp(c);
            int mf, of, df, mr, orr, dr;
            place(r, mf, of, df);
            place(rcbuf, mr, orr, dr);
            const bool use_rc = (mr * of > mf * orr) ||
                                (mr * of == mf * orr && mr > mf);
            const std::string& q = use_rc ? rcbuf : r;
            const int m = use_rc ? mr : mf, ov = use_rc ? orr : of, d = use_rc ? dr : df;
            if (m < 0 || m < static_cast<int>(min_ident * ov)) continue;
            int lo = std::max(0, d), hi = std::min(L, d + static_cast<int>(q.size()));
            for (int p = lo; p < hi; ++p) {
                int b = bidx(q[p - d]);
                if (b >= 0) ++cnt[p][b];
            }
        }
        for (int i = 0; i < L; ++i) {
            int b = bidx(contig[i]);
            int bb = (b >= 0) ? b : 0;
            for (int x = 0; x < 4; ++x) if (cnt[i][x] > cnt[i][bb]) bb = x;
            if (cnt[i][bb] > 0) contig[i] = BASES[bb];
        }
    }
}

// Reads-mode locus clustering (always on, sequence-only): two reads are linked
// into the same locus cluster when their weighted shared k-mer support reaches
// kClusterMinScore. Each shared element k-mer contributes 1/(number of reads
// carrying it), so cross-locus repeat k-mers (present in many reads) barely
// count while same-locus k-mers count strongly. Clusters are assembled
// separately and the best-supported one is kept, so repeat reads cannot
// pollute the contig.
static constexpr double kClusterMinScore = 2.0;

static std::vector<int> cluster_reads_by_shared_kmers(
    const std::vector<std::string>& reads, int k,
    double min_score = kClusterMinScore, bool count_mode = false)
{
    const int n = static_cast<int>(reads.size());
    std::vector<int> parent(n);
    for (int i = 0; i < n; ++i) parent[i] = i;
    if (n < 2) return parent;

    auto find_root = [&](int x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    };

    std::unordered_map<uint32_t, std::vector<int>> occ;
    for (int i = 0; i < n; ++i) {
        const std::string& s = reads[i];
        const int last = static_cast<int>(s.size()) - k;
        for (int p = 0; p <= last; ++p) {
            uint32_t code;
            if (!encode_mer_at(s, static_cast<size_t>(p), static_cast<size_t>(k), code))
                continue;
            auto& lst = occ[code];
            if (!lst.empty() && lst.back() == i) continue;
            lst.push_back(i);
        }
    }

    std::unordered_map<uint64_t, double> pair_score;
    for (const auto& kv : occ) {
        const auto& lst = kv.second;
        if (lst.size() < 2 || lst.size() > 128) continue;
        const double w = count_mode ? 1.0 : 1.0 / static_cast<double>(lst.size());
        for (size_t a = 0; a + 1 < lst.size(); ++a) {
            for (size_t b = a + 1; b < lst.size(); ++b) {
                int x = lst[a], y = lst[b];
                if (x > y) std::swap(x, y);
                pair_score[(static_cast<uint64_t>(x) << 32) |
                           static_cast<uint32_t>(y)] += w;
            }
        }
    }
    for (const auto& kv : pair_score) {
        if (kv.second < min_score) continue;
        int a = find_root(static_cast<int>(kv.first >> 32));
        int b = find_root(static_cast<int>(kv.first & 0xffffffffu));
        if (a != b) parent[b] = a;
    }
    for (int i = 0; i < n; ++i) parent[i] = find_root(i);
    return parent;
}

std::unordered_map<int, std::string> load_cne_refs(const std::string& path) {
    std::unordered_map<int, std::string> refs;
    BufferedLineReader r(path);
    std::string line, seq;
    int id = 0;
    while (r.readline(line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!seq.empty() && id > 0) refs[id] = std::move(seq);
            seq.clear();
            id = std::stoi(line.substr(1));
        } else {
            seq.append(line);
        }
    }
    if (!seq.empty() && id > 0) refs[id] = std::move(seq);
    std::cerr << "  Loaded " << refs.size() << " cne refs\n";
    return refs;
}

// Local SW identity with the evaluation scoring (match +2, mismatch -1,
// linear gap 1/base). Returns matches / alignment length of the best local
// alignment, 0.0 when no positive-scoring alignment exists.
static double sw_identity_2m1(const std::string& ref, const std::string& seq) {
    const int m = static_cast<int>(ref.size());
    const int n = static_cast<int>(seq.size());
    if (m == 0 || n == 0) return 0.0;
    const int match_score = 2, mismatch_score = -1, gap = 1;
    std::vector<std::vector<int>> M(m + 1, std::vector<int>(n + 1, 0));
    std::vector<std::vector<int>> X(m + 1, std::vector<int>(n + 1, 0));
    std::vector<std::vector<int>> Y(m + 1, std::vector<int>(n + 1, 0));
    int best = 0, max_i = 0, max_j = 0, best_state = 0;
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            int d = (ref[i - 1] == seq[j - 1]) ? match_score : mismatch_score;
            X[i][j] = std::max({M[i - 1][j] - gap, X[i - 1][j] - gap, 0});
            Y[i][j] = std::max({M[i][j - 1] - gap, Y[i][j - 1] - gap, 0});
            M[i][j] = std::max({0, M[i - 1][j - 1] + d, X[i][j], Y[i][j]});
            if (M[i][j] > best) { best = M[i][j]; max_i = i; max_j = j; best_state = 0; }
            if (X[i][j] > best) { best = X[i][j]; max_i = i; max_j = j; best_state = 1; }
            if (Y[i][j] > best) { best = Y[i][j]; max_i = i; max_j = j; best_state = 2; }
        }
    }
    if (best == 0) return 0.0;
    int i = max_i, j = max_j, state = best_state;
    int matches = 0, alen = 0;
    while (i > 0 && j > 0) {
        if (state == 0) {
            if (M[i][j] == 0) break;
            int d = (ref[i - 1] == seq[j - 1]) ? match_score : mismatch_score;
            if (M[i][j] == M[i - 1][j - 1] + d) {
                if (ref[i - 1] == seq[j - 1]) ++matches;
                ++alen;
                --i; --j;
                continue;
            }
            if (M[i][j] == X[i][j]) { state = 1; continue; }
            if (M[i][j] == Y[i][j]) { state = 2; continue; }
            break;
        } else if (state == 1) {
            ++alen;
            if (X[i][j] == M[i - 1][j] - gap) state = 0;
            --i;
        } else {
            ++alen;
            if (Y[i][j] == M[i][j - 1] - gap) state = 0;
            --j;
        }
    }
    return alen > 0 ? static_cast<double>(matches) / alen : 0.0;
}

struct GenomeResult {
    int ele_id;
    std::string display_id;
    std::string seq;
};

static void process_genome_element(
    const std::vector<Read>& all_reads, size_t begin, size_t end,
    const Args& args, const MerQueryManager& mqm,
    const PanelScorer* psc,
    const std::unordered_map<int, std::string>* cne_refs,
    std::vector<GenomeResult>& out)
{
    const int cur_ele = all_reads[begin].ele_id;
    struct EleRead {
        Read read;
        std::string seq;
        int seen = 0;   // distinct element mers in the read
        int pos = 0;    // matched element-mer positions in the read
    };
    std::vector<EleRead> ele_reads;
    ele_reads.reserve(end - begin);
    size_t idx = begin;
    while (idx < end) {
        const auto& r = all_reads[idx];
        std::string local_seq = r.seq;
        normalize_read_seq(local_seq, r.strand);
        // Repeat ratio check for genome mode; the distinct element-mer count
        // doubles as a per-read score bound for branch-and-bound below.
        {
            thread_local tsl::robin_set<uint32_t> seen;
            seen.clear();
            int total = 0;
            int ms = mqm.get_mer_size();
            for (int k = 0; k + ms <= (int)local_seq.size(); ++k) {
                uint32_t code;
                if (encode_mer_at(local_seq, k, ms, code)) {
                    auto it = mqm.compressed_mer_query.find(code);
                    if (it != mqm.compressed_mer_query.end() && mqm.has_id(it->second, cur_ele)) {
                        ++total;
                        seen.insert(code);
                    }
                }
            }
            double ratio = (double)total / std::max(1, (int)seen.size());
            if (ratio > args.repeat_ratio) {
                ++idx;
                continue;
            }
            ele_reads.push_back(EleRead{r, std::move(local_seq),
                                        static_cast<int>(seen.size()), total});
        }
        ++idx;
    }
    if (ele_reads.empty()) return;

    // Cluster by locus
    std::sort(ele_reads.begin(), ele_reads.end(),
        [](const EleRead& a, const EleRead& b) {
            if (a.read.chr != b.read.chr) return a.read.chr < b.read.chr;
            if (a.read.strand != b.read.strand) return a.read.strand > b.read.strand;
            return a.read.start < b.read.start;
        });

    struct GenomeCand {
        std::string seq;
        std::string did;
        ContigScore sc{0, 0};
        size_t support = 0;
        double pid = -1.0;
        size_t pos = 0;
    };
    std::vector<GenomeCand> cands;
    cands.reserve(16);
    double best_score = -1.0;
    size_t best_pos = 0;
    const bool panels_on = psc != nullptr && psc->enabled();

    auto consider = [&](std::string contig, int cluster_start,
                        const EleRead& first_r, const EleRead& last_r,
                        size_t pos, size_t support) {
        if (contig.empty()) return;
        if (args.trim) trim_contig(contig, mqm, cur_ele, args.trim_margin);
        if (contig.size() < 20) return;
        ContigScore sc = score_contig(contig, mqm, cur_ele);
        bool repeat_chimera = sc.hits >= 10 && sc.distinct > 0 &&
            static_cast<double>(sc.hits) / sc.distinct > 2.0;
        if (repeat_chimera) return;
        char s = first_r.read.strand == 1 ? '+' : '-';
        std::string did = std::to_string(cur_ele) + "." +
            first_r.read.chr + ":" +
            std::to_string(cluster_start) + "-" +
            std::to_string(last_r.read.end) + "(" + s + ")";
        GenomeCand cand;
        cand.seq = std::move(contig);
        cand.did = std::move(did);
        cand.sc = sc;
        cand.support = support;
        cand.pos = pos;
        if (panels_on) cand.pid = psc->score(cur_ele, cand.seq, 0.5).max_id;
        cands.push_back(std::move(cand));
        if (sc.hits > best_score || (sc.hits == best_score && pos < best_pos)) {
            best_score = static_cast<double>(sc.hits);
            best_pos = pos;
        }
    };

    // Sample one cluster's reads down to the per-locus budget.
    auto sample_cluster = [&](const std::vector<EleRead>& er, size_t cb, size_t ce) {
        size_t n_reads = ce - cb;
        std::vector<std::string> seqs;
        seqs.reserve(std::min<size_t>(n_reads, args.max_reads));
        if ((int)n_reads <= args.max_reads) {
            for (size_t k = cb; k < ce; ++k) seqs.push_back(er[k].seq);
            return seqs;
        }
        int n = args.max_reads;
        for (int i = 0; i < n; ++i) {
            size_t j = (n == 1)
                ? 0
                : static_cast<size_t>((double)i * (n_reads - 1) / (n - 1) + 0.5);
            seqs.push_back(er[cb + j].seq);
        }
        return seqs;
    };

    // Cluster ranges with a cheap, valid score upper bound. A contig's hits
    // are bounded by its length (span - 12); the repeat guard also rejects any
    // contig with hits >= 10 and hits/distinct > 2, so accepted contigs have
    // hits <= max(9, 2 * distinct), and distinct mers cannot exceed the sum of
    // per-read distinct element mers in the cluster.
    struct ClusterRange {
        size_t begin, end;
        int ub;
        int span;
        size_t pos;   // positional (legacy) order, used for tie-breaks
    };
    std::vector<ClusterRange> ranges;
    size_t cb = 0;
    for (size_t j = 1; j <= ele_reads.size(); ++j) {
        bool cut = (j == ele_reads.size());
        if (!cut) {
            int gap = ele_reads[j].read.start - ele_reads[j - 1].read.end;
            if (gap > args.max_loci_gap ||
                ele_reads[j].read.chr != ele_reads[j - 1].read.chr ||
                ele_reads[j].read.strand != ele_reads[j - 1].read.strand) {
                cut = true;
            }
        }
        if (cut) {
            long long sum_seen = 0, sum_pos = 0;
            for (size_t k = cb; k < j; ++k) {
                sum_seen += ele_reads[k].seen;
                sum_pos += ele_reads[k].pos;
            }
            int span = ele_reads[j - 1].read.end - ele_reads[cb].read.start;
            long long guard_ub = std::max<long long>(9, 2 * sum_seen);
            long long ub = std::min<long long>(
                guard_ub, std::min<long long>(sum_pos, std::max(0, span - 12)));
            ranges.push_back(ClusterRange{cb, j, static_cast<int>(ub), span, ranges.size()});
            cb = j;
        }
    }

    // Branch and bound: assemble by descending upper bound and stop once the
    // bound can no longer match the best score. Clusters whose bound equals the
    // current best can only tie it, so they are deferred and later examined in
    // legacy positional order (on equal hits the earliest locus wins, exactly
    // like the exhaustive positional scan).
    std::sort(ranges.begin(), ranges.end(),
              [](const ClusterRange& a, const ClusterRange& b) { return a.ub > b.ub; });
    bool debug_ele = args.debug_eles.count(cur_ele) > 0;
    auto try_range = [&](const ClusterRange& cr) {
        const size_t n_reads = cr.end - cr.begin;
        if (n_reads < 2) {
            // A de Bruijn path on a single window is the window itself, and
            // edge filtering cannot remove spurious/paralogous windows: with
            // min-count=1 every lone window would become a claim. Mirror the
            // reads-mode single-sequence fallback instead: keep the window
            // only when it carries enough distinct element mers itself;
            // weaker ones are left to the T3 rescue handoff.
            const EleRead& er = ele_reads[cr.begin];
            if (er.seen < kMinFallbackDistinct) return;
            consider(er.seq, er.read.start, er, er, cr.pos, 1);
            return;
        }
        std::vector<std::string> reads_for_graph = sample_cluster(ele_reads, cr.begin, cr.end);
        std::string contig = debruijn_assemble(reads_for_graph, args.kmer, args.min_count,
                                               mqm, cur_ele, nullptr, nullptr, debug_ele);
        consider(std::move(contig), ele_reads[cr.begin].read.start,
                 ele_reads[cr.begin], ele_reads[cr.end - 1], cr.pos,
                 cr.end - cr.begin);
    };

    std::vector<const ClusterRange*> ties;
    for (const auto& cr : ranges) {
        if (static_cast<double>(cr.ub) < best_score) break;
        if (cr.span < 20) continue;
        if (static_cast<double>(cr.ub) == best_score) {
            ties.push_back(&cr);
            continue;
        }
        try_range(cr);
    }
    std::sort(ties.begin(), ties.end(),
              [](const ClusterRange* a, const ClusterRange* b) { return a->pos < b->pos; });
    for (const ClusterRange* cr : ties) {
        if (cr->pos >= best_pos) break;
        if (static_cast<double>(cr->ub) < best_score) continue;
        try_range(*cr);
    }

    if (!cands.empty()) {
        size_t best_i = 0;
        for (size_t i = 1; i < cands.size(); ++i) {
            if (cands[i].sc.hits > cands[best_i].sc.hits ||
                (cands[i].sc.hits == cands[best_i].sc.hits &&
                 cands[i].pos < cands[best_i].pos)) {
                best_i = i;
            }
        }
        GenomeCand& bc = cands[best_i];
        if (panels_on && bc.sc.distinct < args.min_output_distinct) {
            bool keep = false;
            if (panels_on && bc.sc.distinct >= args.rescue_min_distinct &&
                target_span(bc.seq, mqm, cur_ele) >= args.rescue_min_span) {
                double pid = bc.pid;
                if (pid < 0.0) pid = psc->identity(cur_ele, bc.seq, 0.5);
                if (pid >= args.rescue_tau) keep = true;
            }
            if (!keep && bc.sc.distinct >= kMinOutputDistinctLow &&
                bc.support >= kMinOutputSupportLow && cne_refs != nullptr &&
                !cne_refs->empty()) {
                auto it = cne_refs->find(cur_ele);
                if (it != cne_refs->end() &&
                    sw_identity_2m1(it->second, bc.seq) >= min_cne_ref_identity()) {
                    keep = true;
                }
            }
            if (!keep) bc.seq.clear();
        }
        if (!bc.seq.empty())
            out.push_back(GenomeResult{cur_ele, std::move(bc.did), std::move(bc.seq)});
    } else if (rescue_handoff_enabled()) {
        // T3 rescue handoff: validate already SW-verified these windows against
        // the element panel; mer-based gates do not apply to them. Emit the
        // best-scoring window as the element's claim.
        const Read* best_r = nullptr;
        double best_rs = 0.0;
        for (size_t i = begin; i < end; ++i) {
            const auto& r = all_reads[i];
            if (r.rescue > best_rs) {
                best_rs = r.rescue;
                best_r = &r;
            }
        }
        if (best_r != nullptr) {
            std::string seq = best_r->seq;
            normalize_read_seq(seq, best_r->strand);
            if (seq.size() >= 20) {
                char s = best_r->strand == 1 ? '+' : '-';
                std::string did = std::to_string(cur_ele) + "." + best_r->chr + ":" +
                    std::to_string(best_r->start) + "-" + std::to_string(best_r->end) +
                    "(" + s + ")";
                out.push_back(GenomeResult{cur_ele, std::move(did), std::move(seq)});
            }
        }
    }
}

int main(int argc, char* argv[]) {
    try {
        Args args = parse_args(argc, argv);

        if (args.help_requested || args.input_dir.empty()) {
            std::cerr << "Usage: " << argv[0]
                      << " <input_dir> --mers <mers_table> [options]\n"
                      << "\nOptions:\n"
                      << "  --mers <file>         Mers table (required)\n"
                      << "  -o, --output <file>   Output FASTA (default: assembled.fasta)\n"
                      << "  -k, --kmer <n>        K-mer length, single-k mode (default: 19)\n"
                      << "  --multi-k <k1,k2,..>  Reads mode: assemble each element with every k\n"
                      << "                        and keep the best-scoring contig\n"
                      << "                        (default: 27; -k or --no-multi-k disables)\n"
                      << "  --no-multi-k          Force single-k mode\n"
                      << "  --min-count <n>       Min (k+1)-mer occurrence in graph (default: 1)\n"
                      << "  --max-reads <n>       Max reads per element/locus (default: 200)\n"
                      << "  --max-loci-gap <n>    Max gap for locus clustering (default: 5000)\n"
                      << "  --max-assign <n>      Load-time cap on (element, locus) assignments per\n"
                      << "                        k-mer, must match validate (default: 3; 0 = keep all)\n"
                      << "  -t, --threads <n>     Worker threads for genome-mode assembly (default: 4)\n"
                      << "  --trim                Trim contigs to confident k-mer region (default: on)\n"
                      << "  --no-trim             Disable trimming\n"
                      << "  --trim-margin <n>     bp kept beyond first/last element k-mer (default: 10)\n"
                      << "  --no-polish           Disable reads-mode post-assembly polish\n"
                      << "  --repeat-ratio <f>    Max total/unique mers per read (default: 1.2)\n"
                       << "  --snp                 Scan for candidate SNPs from validated reads\n"
                       << "  --cne-ref <file>      CNE reference FASTA (ref.fa): identity verification of\n"
                      << "                        low-distinct claims (reads mode)\n";
            return 1;
        }

        std::cerr << "Loading mers table from " << args.mers_file << " ...\n";
        MerQueryManager mqm;
        mqm.load_from_file(args.mers_file, args.max_assign);
        int mer_size = mqm.get_mer_size();
        std::cerr << "  Loaded " << mqm.size() << " mers, k=" << mer_size << "\n";

        // Load CNE refs for SW-based trimming
        std::unordered_map<int, std::string> cne_refs;
        if (!args.cne_ref_path.empty()) {
            std::cerr << "Loading CNE refs from " << args.cne_ref_path << " ...\n";
            cne_refs = load_cne_refs(args.cne_ref_path);
        }

        // Optional leave-one-out panels for cluster selection / low-support
        // rescue. The strict path is untouched when this is empty.
        std::unordered_map<int, std::vector<std::string>> panels;
        parasail_matrix_t* panel_mat = nullptr;
        PanelScorer psc;
        if (args.rescue_enabled && !args.panels_path.empty()) {
            std::cerr << "Loading rescue panels from " << args.panels_path << " ...\n";
            panel_core::load_panels(args.panels_path, 3, panels);
            size_t nseq = 0;
            for (const auto& kv : panels) nseq += kv.second.size();
            panel_mat = parasail_matrix_create("ACGTN", 2, -1);
            psc.panels = &panels;
            psc.mat = panel_mat;
            std::cerr << "  panels: " << panels.size() << " elements, " << nseq
                      << " sequences; rescue tau=" << args.rescue_tau << "\n";
        }

        // Phase 1: Read all Assemble.*.reads files
        std::vector<Read> all_reads;
        bool is_genome = false;

        if (!fs::exists(args.input_dir)) {
            throw std::runtime_error("Input directory not found: " + args.input_dir);
        }

        for (const auto& entry : fs::directory_iterator(args.input_dir)) {
            auto path = entry.path().string();
            if (path.find("Assemble.") == std::string::npos) continue;
            if (path.rfind(".reads") != path.size() - 6) continue;

            BufferedLineReader reader(path);
            std::string line;
            while (reader.readline(line)) {
                if (line.empty()) continue;
                // Parse 4 columns
                size_t t1 = line.find('\t');
                if (t1 == std::string::npos) continue;
                size_t t2 = line.find('\t', t1 + 1);
                if (t2 == std::string::npos) continue;
                size_t t3 = line.find('\t', t2 + 1);
                if (t3 == std::string::npos) continue;

                std::string seq_id = line.substr(0, t1);
                int strand = std::stoi(line.substr(t1 + 1, t2 - t1 - 1));
                int ele_id = std::stoi(line.substr(t2 + 1, t3 - t2 - 1));
                // Optional 5th column: rescue score from validate (tab-separated).
                size_t t4 = line.find('\t', t3 + 1);
                std::string seq = (t4 == std::string::npos)
                    ? line.substr(t3 + 1)
                    : line.substr(t3 + 1, t4 - t3 - 1);

                if (seq.empty()) continue;

                Read read;
                read.seq_id = seq_id;
                read.strand = strand;
                read.ele_id = ele_id;
                read.seq = seq;
                if (t4 != std::string::npos) {
                    try { read.rescue = std::stod(line.substr(t4 + 1)); }
                    catch (...) { read.rescue = 0.0; }
                }

                GenomePos gp;
                if (!is_genome && parse_genome_id(seq_id, gp)) {
                    is_genome = true;
                }
                if (is_genome) {
                    if (parse_genome_id(seq_id, gp)) {
                        read.chr = gp.chr;
                        read.start = gp.start;
                        read.end = gp.end;
                    }
                }

                all_reads.push_back(std::move(read));
            }
        }

        // Drop low-scoring T3 rescue reads unless they reach the configured
        // score threshold.
        if (args.min_rescue_score > 0.0 && !is_genome) {
            const double thr = args.min_rescue_score;
            size_t before = all_reads.size();
            all_reads.erase(std::remove_if(all_reads.begin(), all_reads.end(),
                [&](const Read& r) {
                    return r.rescue > 0.0 && r.rescue < thr;
                }), all_reads.end());
            std::cerr << "  rescue filter: kept " << all_reads.size() << "/" << before
                      << " reads (score>=" << thr << ")\n";
        }
        std::cerr << "  Loaded " << all_reads.size() << " reads\n";
        if (all_reads.empty()) {
            std::cerr << "No reads found, output empty.\n";
            std::ofstream ofs(args.output);
            return 0;
        }

        // Phase 2: Sort
        if (is_genome) {
            std::sort(all_reads.begin(), all_reads.end(),
                [](const Read& a, const Read& b) {
                    if (a.ele_id != b.ele_id) return a.ele_id < b.ele_id;
                    if (a.chr != b.chr) return a.chr < b.chr;
                    return a.start < b.start;
                });
        } else {
            std::sort(all_reads.begin(), all_reads.end(),
                [](const Read& a, const Read& b) {
                    if (a.ele_id != b.ele_id) return a.ele_id < b.ele_id;
                    return a.seq_id < b.seq_id;
                });
        }

        // Phase 3: Group and assemble
        std::vector<std::pair<std::string, std::string>> results;

        // Open SNP output if requested (reads mode only)
        std::ofstream snp_stream, gfa_stream;
        std::ofstream* snp_out_ptr = nullptr;
        std::ofstream* gfa_out_ptr = nullptr;
        if (args.output_snp && !is_genome) {
            snp_stream.open("variants.tsv");
            snp_stream << "ele_id\tpos\tref\talt\tref_cov\talt_cov\talt_freq\tbranch_len\ttype\n";
            snp_out_ptr = &snp_stream;
            gfa_stream.open("snp_elements.gfa");
            gfa_stream << "H\tVN:Z:1.0\n";
            gfa_out_ptr = &gfa_stream;
        }

        // Per-element reads-mode worker: returns the claim sequence ("" = none).
        // Pure w.r.t. shared state (mqm/cne_refs are read-only), so it can run
        // from several threads; the only shared resource is the optional SNP/GFA
        // stream pair, which is why the SNP path stays sequential.
        auto process_reads_element = [&](const std::vector<Read>& pool,
                                         size_t begin, size_t end,
                                         std::ofstream* snp_out_ptr,
                                         std::ofstream* gfa_out_ptr) -> std::string {
                size_t idx = begin;
                int cur_ele = pool[idx].ele_id;

                // Reads mode: one group per element
                std::vector<std::string> reads;
                reads.reserve(args.max_reads);
                std::vector<double> rescues;   // parallel to reads: T3 rescue score
                rescues.reserve(args.max_reads);
                int cnt = 0;
                while (idx < end && cnt < args.max_reads) {
                    const auto& r = pool[idx];
                    std::string local_seq = r.seq;
                    normalize_read_seq(local_seq, r.strand);
                    // Repeat ratio check
                    {
                        thread_local tsl::robin_set<uint32_t> seen;
                        seen.clear();
                        int total = 0;
                        int ms = mqm.get_mer_size();
                        for (int k = 0; k + ms <= (int)local_seq.size(); ++k) {
                            uint32_t code;
                            if (encode_mer_at(local_seq, k, ms, code)) {
                                auto it = mqm.compressed_mer_query.find(code);
                                if (it != mqm.compressed_mer_query.end() && mqm.has_id(it->second, cur_ele)) {
                                    ++total;
                                    seen.insert(code);
                                }
                            }
                        }
                        double ratio = (double)total / std::max(1, (int)seen.size());
                        if (ratio > args.repeat_ratio) {
                            ++idx;
                            continue;
                        }
                    }
                    reads.push_back(std::move(local_seq));
                    rescues.push_back(r.rescue);
                    ++idx;
                    ++cnt;
                }

                // Best-cluster assembly: split the element's reads into
                // locus-pure clusters, assemble each, trim, apply the
                // repeat-chimera guard and retain the candidates. When panels
                // are available, choose the panel-best candidate among those
                // with competitive hit counts (a paralog often wins on raw
                // element-mer hits but loses on ortholog identity).
                struct Cand {
                    std::string seq;
                    ContigScore sc{0, 0};
                    size_t support = 0;
                    double pid = -1.0;
                };
                std::vector<Cand> cands;
                cands.reserve(16);
                std::string contig;
                ContigScore contig_sc{0, 0};
                double contig_pid = -1.0;
                size_t best_support = 0;
                int last_second_hits = -1;
                int last_total_reads = 0;
                int last_groups = 0;
                auto dbg_bail = [&](const char* tag, long v) {
                    if (asm_debug_enabled())
                        std::fprintf(stderr, "ASMDBG2\t%d\t%s\t%ld\t%zu\n",
                                     cur_ele, tag, v, reads.size());
                };
                int dbg_ng = 0, dbg_cand = 0, dbg_trim20 = 0, dbg_chim = 0;
                bool dbg_floor_cleared = false;
                if (reads.size() >= 2) {
                    std::vector<int> ks = args.multi_k.empty()
                        ? std::vector<int>{args.kmer} : args.multi_k;
                    int total_groups = 0;
                    for (int kk : ks) {
                        std::vector<int> labels =
                            cluster_reads_by_shared_kmers(reads, kk);
                        std::unordered_map<int, std::vector<int>> groups;
                        for (int i = 0; i < static_cast<int>(labels.size()); ++i)
                            groups[labels[i]].push_back(i);
                        total_groups += static_cast<int>(groups.size());
                        for (const auto& kv : groups) {
                            if (kv.second.size() < 2) continue;
                            ++dbg_ng;
                            std::vector<std::string> sub;
                            sub.reserve(kv.second.size());
                            for (int i : kv.second) sub.push_back(reads[i]);
                            Cand cd;
                            std::string c = debruijn_assemble(sub, kk, args.min_count,
                                                              mqm, cur_ele, snp_out_ptr,
                                                              gfa_out_ptr, false);
                            if (c.empty()) continue;
                            ++dbg_cand;
                            if (args.trim) trim_contig(c, mqm, cur_ele, args.trim_margin);
                            if (c.size() < 20) continue;
                            ++dbg_trim20;
                            ContigScore sc = score_contig(c, mqm, cur_ele);
                            bool repeat_chimera = sc.hits >= 10 && sc.distinct > 0 &&
                                static_cast<double>(sc.hits) / sc.distinct > 2.0;
                            if (repeat_chimera) continue;
                            ++dbg_chim;
                            cd.seq = std::move(c);
                            cd.sc = sc;
                            cd.support = kv.second.size();
                            cands.push_back(std::move(cd));
                        }
                    }
                    last_total_reads = static_cast<int>(reads.size());
                    last_groups = total_groups;
                }

                if (!cands.empty()) {
                    size_t best_i = 0;
                    for (size_t i = 1; i < cands.size(); ++i) {
                        const Cand& a = cands[i];
                        const Cand& b = cands[best_i];
                        auto path_key = [&](const Cand& c) -> double {
                            return c.sc.hits;
                        };
                        double ka = path_key(a), kb = path_key(b);
                        if (ka > kb ||
                            (ka == kb && a.support > b.support) ||
                            (ka == kb && a.support == b.support && a.seq.size() > b.seq.size())) {
                            best_i = i;
                        }
                    }
                    size_t sel = best_i;
                    if (psc.enabled()) {
                        for (auto& cd : cands) {
                            PanelScore ps = psc.score(cur_ele, cd.seq, 0.5);
                            cd.pid = ps.median_id;   // robust across panel species
                        }
                        // Conservative rescue: only override a hits-best
                        // candidate when that candidate itself fails the panel
                        // threshold. This avoids replacing already-correct
                        // high-identity clusters (the old unconditional
                        // panel-best rule lost true positives).
                        if (cands[best_i].pid < args.rescue_tau) {
                            const double hits_max = static_cast<double>(cands[best_i].sc.hits);
                            size_t best_pid = best_i;
                            for (size_t i = 0; i < cands.size(); ++i) {
                                if (static_cast<double>(cands[i].sc.hits) < 0.5 * hits_max) continue;
                                if (cands[i].pid > cands[best_pid].pid + 1e-9) best_pid = i;
                            }
                            if (best_pid != best_i && cands[best_pid].pid >= args.rescue_tau)
                                sel = best_pid;
                        }
                    }
                    contig = std::move(cands[sel].seq);
                    contig_sc = cands[sel].sc;
                    contig_pid = cands[sel].pid;
                    best_support = cands[sel].support;
                    int sh = -1;
                    for (size_t i = 0; i < cands.size(); ++i)
                        if (i != sel && static_cast<int>(cands[i].sc.hits) > sh)
                            sh = static_cast<int>(cands[i].sc.hits);
                    last_second_hits = sh;
                }
                // Output evidence floor: a claim needs at least this many
                // distinct element mers, otherwise it is a short/repeat junk
                // fragment (precision parity with genome mode). A contig that
                // fails the floor counts as "no acceptable output".
                // `final_sc` is the score of the sequence that is actually
                // emitted. The cluster winner was already scored at selection
                // time and its sequence is unchanged here, so reuse it; the
                // fallback/rescue paths change the sequence and rescore below.
                ContigScore final_sc{0, 0};
                bool final_sc_valid = false;
                if (!contig.empty()) {
                    final_sc = contig_sc;
                    final_sc_valid = true;
                    if (final_sc.distinct < args.min_output_distinct) {
                        bool keep = false;
                        if (psc.enabled() && final_sc.distinct >= args.rescue_min_distinct &&
                            target_span(contig, mqm, cur_ele) >= args.rescue_min_span) {
                            double pid = contig_pid;
                            if (pid < 0.0) pid = psc.identity(cur_ele, contig, 0.5);
                            if (pid >= args.rescue_tau) keep = true;
                        }
                        if (!keep && final_sc.distinct >= kMinOutputDistinctLow &&
                            best_support >= kMinOutputSupportLow && !cne_refs.empty()) {
                            auto it = cne_refs.find(cur_ele);
                            if (it != cne_refs.end() &&
                                sw_identity_2m1(it->second, contig) >= min_cne_ref_identity()) {
                                keep = true;
                            }
                        }
                        if (!keep) {
                            contig.clear();
                            final_sc_valid = false;
                            dbg_floor_cleared = true;
                        }
                    }
                }
                if (args.rescue_lowconf && final_sc_valid && !reads.empty()) {
                    const int cur_dist = final_sc.distinct;
                    int best_i = -1, best_dist = -1;
                    std::vector<int> labels = cluster_reads_by_shared_kmers(reads, mqm.get_mer_size());
                    std::unordered_map<int, std::vector<int>> groups;
                    for (int i = 0; i < static_cast<int>(labels.size()); ++i)
                        groups[labels[i]].push_back(i);
                    for (const auto& kv : groups) {
                        if (static_cast<int>(kv.second.size()) < args.rescue_min_support) continue;
                        for (int i : kv.second) {
                            ContigScore sc = score_contig(reads[i], mqm, cur_ele);
                            if (sc.distinct > best_dist) { best_dist = sc.distinct; best_i = i; }
                        }
                    }
                    if (best_i >= 0 && best_dist > cur_dist) {
                        contig = reads[best_i];
                        final_sc = score_contig(contig, mqm, cur_ele);
                        final_sc_valid = true;
                        dbg_floor_cleared = false;
                    }
                }
                if (asm_debug_enabled()) {
                    std::fprintf(stderr,
                                 "ASMDBG2\t%d\tnr=%zu\tng=%d\tcand=%d\ttrim20=%d\tchim=%d\tcontig=%zu\tsup=%zu\tfloor=%d\tsh=%d\n",
                                 cur_ele, reads.size(), dbg_ng, dbg_cand, dbg_trim20,
                                 dbg_chim, contig.size(), best_support,
                                 dbg_floor_cleared ? 1 : 0, last_second_hits);
                }

                bool rescued = false;
                if (contig.empty()) {
                    // Fallback: the element's verified reads are still
                    // evidence, so emit the single read carrying the most
                    // distinct element mers (a cheap, alignment-free fidelity
                    // proxy: a true-locus read covers many different element
                    // mers, repeat-derived noise covers few). Same semantics as
                    // genome mode: no claim unless the floor is cleared.
                    int best_dist = -1;
                    size_t best_idx = 0;
                    for (size_t i = 0; i < reads.size(); ++i) {
                        ContigScore sc = score_contig(reads[i], mqm, cur_ele);
                        if (sc.distinct > best_dist) {
                            best_dist = sc.distinct;
                            best_idx = i;
                        }
                    }
                    if (best_dist >= args.min_fallback_distinct) {
                        contig = std::move(reads[best_idx]);
                        if (args.trim) trim_contig(contig, mqm, cur_ele, args.trim_margin);
                        if (contig.size() < 20) { dbg_bail("fb_short", best_dist); return std::string(); }
                    } else if (psc.enabled()) {
                        // Panel rescue: a single verified read can be the only
                        // evidence for a short/low-complexity element. Require
                        // a strict panel identity and a minimal mer footprint;
                        // emit the read untrimmed (trimming a sparse-mer read
                        // would cut it to nothing).
                        size_t ri = reads.size();
                        double best_pid = -1.0;
                        for (size_t i = 0; i < reads.size(); ++i) {
                            ContigScore sc = score_contig(reads[i], mqm, cur_ele);
                            if (sc.distinct < args.rescue_min_distinct) continue;
                            if (target_span(reads[i], mqm, cur_ele) < args.rescue_min_span) continue;
                            double pid = psc.identity(cur_ele, reads[i], 0.5);
                            if (pid > best_pid) { best_pid = pid; ri = i; }
                        }
                        if (asm_debug_enabled())
                            std::fprintf(stderr, "ASMDBG2\t%d\tfallback\tbest_dist=%d\tbest_pid=%.4f\tnr=%zu\n",
                                         cur_ele, best_dist, best_pid, reads.size());
                        if (ri < reads.size() && best_pid >= args.rescue_tau) {
                            contig = std::move(reads[ri]);
                            rescued = true;
                        } else {
                            // T3 rescue handoff: validate already SW-verified
                            // these reads against the element panel.
                            if (!rescue_handoff_enabled()) { dbg_bail("no_handoff", best_dist); return std::string(); }
                            size_t rj = reads.size();
                            double best_rs = 0.0;
                            for (size_t i = 0; i < reads.size(); ++i) {
                                if (rescues[i] > best_rs) {
                                    best_rs = rescues[i];
                                    rj = i;
                                }
                            }
                            if (rj == reads.size()) { dbg_bail("no_rescue", best_dist); return std::string(); }
                            contig = std::move(reads[rj]);
                            rescued = true;
                            if (contig.size() < 20) { dbg_bail("rc_short", (long)best_rs); return std::string(); }
                        }
                    } else {
                        // T3 rescue handoff: validate already SW-verified these
                        // reads against the element panel, so the mer-distinct
                        // floor does not apply. Emit the best-scoring one as-is
                        // (mer-based trimming would cut it to nothing).
                        if (!rescue_handoff_enabled()) { dbg_bail("no_handoff", best_dist); return std::string(); }
                        size_t ri = reads.size();
                        double best_rs = 0.0;
                        for (size_t i = 0; i < reads.size(); ++i) {
                            if (rescues[i] > best_rs) {
                                best_rs = rescues[i];
                                ri = i;
                            }
                        }
                        if (ri == reads.size()) { dbg_bail("no_rescue", best_dist); return std::string(); }
                        contig = std::move(reads[ri]);
                        rescued = true;
                        if (contig.size() < 20) { dbg_bail("rc_short", (long)best_rs); return std::string(); }
                    }
                }

                if (!final_sc_valid) {
                    if (rescued) {
                        // Rescue claims bypass the mer-distinct floor; the
                        // score is only needed for the debug dump.
                        if (asm_debug_enabled())
                            final_sc = score_contig(contig, mqm, cur_ele);
                    } else {
                        final_sc = score_contig(contig, mqm, cur_ele);
                        if (final_sc.distinct < args.min_output_distinct) {
                            dbg_bail("out_floor", final_sc.distinct);
                            return std::string();
                        }
                    }
                }
                if (asm_debug_enabled()) {
                    std::fprintf(stderr, "ASMDBG\t%d\t%zu\t%d\t%d\t%zu\t%d\t%d\t%d\n",
                                 cur_ele, contig.size(), final_sc.hits,
                                 final_sc.distinct, best_support,
                                 last_total_reads, last_groups, last_second_hits);
                }
                if (args.polish) polish_contig(contig, reads);
                return contig;
        };

        if (!is_genome) {
            // One range per element so the work can fan out across threads.
            std::vector<std::pair<size_t, size_t>> ranges;
            {
                size_t i = 0;
                while (i < all_reads.size()) {
                    size_t j = i + 1;
                    while (j < all_reads.size() && all_reads[j].ele_id == all_reads[i].ele_id) ++j;
                    ranges.emplace_back(i, j);
                    i = j;
                }
            }
            std::vector<std::pair<int, std::string>> eres;
            eres.reserve(ranges.size());
            const int nthreads = std::max(1, args.threads);
            if (args.output_snp || nthreads == 1 || ranges.size() <= 1) {
                for (const auto& r : ranges) {
                    std::string c = process_reads_element(all_reads, r.first, r.second, snp_out_ptr, gfa_out_ptr);
                    if (!c.empty()) eres.emplace_back(all_reads[r.first].ele_id, std::move(c));
                }
            } else {
                std::atomic<size_t> next_range{0};
                std::vector<std::vector<std::pair<int, std::string>>> per_thread(nthreads);
                std::vector<std::thread> workers;
                workers.reserve(nthreads);
                for (int t = 0; t < nthreads; ++t) {
                    workers.emplace_back([&, t]() {
                        while (true) {
                            size_t ri = next_range.fetch_add(1);
                            if (ri >= ranges.size()) break;
                            const auto rr = ranges[ri];
                            std::string c = process_reads_element(all_reads, rr.first, rr.second, nullptr, nullptr);
                            if (!c.empty())
                                per_thread[t].emplace_back(all_reads[rr.first].ele_id, std::move(c));
                        }
                    });
                }
                for (auto& w : workers) w.join();
                for (auto& v : per_thread)
                    for (auto& p : v) eres.push_back(std::move(p));
            }
            std::sort(eres.begin(), eres.end(),
                [](const std::pair<int, std::string>& a, const std::pair<int, std::string>& b) {
                    return a.first < b.first;
                });
            std::map<int, std::string> final_map;
            for (auto& p : eres) final_map.emplace(p.first, std::move(p.second));

            for (auto& kv : final_map)
                results.emplace_back(std::to_string(kv.first), std::move(kv.second));
        } else {
            // Element blocks: all_reads is sorted by ele_id
            std::vector<std::pair<size_t, size_t>> blocks;
            {
                size_t i = 0;
                while (i < all_reads.size()) {
                    size_t j = i + 1;
                    while (j < all_reads.size() && all_reads[j].ele_id == all_reads[i].ele_id) ++j;
                    blocks.emplace_back(i, j);
                    i = j;
                }
            }
            std::vector<GenomeResult> grs;
            const int nthreads = std::max(1, args.threads);
            if (nthreads == 1) {
                for (const auto& b : blocks) {
                    process_genome_element(all_reads, b.first, b.second, args, mqm, &psc, &cne_refs, grs);
                }
            } else {
                std::atomic<size_t> next_block{0};
                std::vector<std::vector<GenomeResult>> per_thread(nthreads);
                std::vector<std::thread> workers;
                workers.reserve(nthreads);
                for (int t = 0; t < nthreads; ++t) {
                    workers.emplace_back([&, t]() {
                        while (true) {
                            size_t bi = next_block.fetch_add(1);
                            if (bi >= blocks.size()) break;
                            process_genome_element(all_reads, blocks[bi].first, blocks[bi].second,
                                                   args, mqm, &psc, &cne_refs, per_thread[t]);
                        }
                    });
                }
                for (auto& w : workers) w.join();
                for (auto& v : per_thread) {
                    for (auto& r : v) grs.push_back(std::move(r));
                }
            }
            std::sort(grs.begin(), grs.end(),
                [](const GenomeResult& a, const GenomeResult& b) { return a.ele_id < b.ele_id; });
            for (auto& r : grs) {
                results.emplace_back(std::move(r.display_id), std::move(r.seq));
            }
        }

        // Phase 4: Write output
        std::ofstream ofs(args.output);
        if (!ofs) {
            throw std::runtime_error("Cannot open output file: " + args.output);
        }
        for (const auto& [seq_id, seq] : results) {
            ofs << ">" << seq_id << "\n" << seq << "\n";
        }

        std::cerr << "Done. " << results.size() << " sequences written to " << args.output << "\n";
        if (panel_mat) parasail_matrix_free(panel_mat);
        return 0;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
