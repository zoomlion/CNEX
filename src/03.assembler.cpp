#include <hip/debruijn_core.hpp>
#include <hip/validator_core.hpp>

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
#include <regex>
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
    int min_count = 1;
    int max_reads = 200;
    int max_loci_gap = 5000;
    int max_assign = 3;   // load-time cap; 0 = keep all
    int threads = 4;
    tsl::robin_set<int> debug_eles;
    bool trim = true;
    bool trim_set = false;
    bool no_trim_set = false;
    double repeat_ratio = 1.2;
    bool output_snp = false;
    bool help_requested = false;
    std::string cne_ref_path;
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
        } else if (arg == "--repeat-ratio") {
            if (++i >= argc) throw std::runtime_error("Missing value for --repeat-ratio");
            args.repeat_ratio = std::stod(argv[i]);
        } else if (arg == "--snp") {
            args.output_snp = true;
        } else if (arg == "--cne-ref") {
            if (++i >= argc) throw std::runtime_error("Missing value for --cne-ref");
            args.cne_ref_path = argv[i];
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

const std::regex GENOME_ID_PATTERN(R"(^(\S+):(\d+)-(\d+)$)");

struct GenomePos {
    std::string chr;
    int start;
    int end;
};

bool parse_genome_id(const std::string& seq_id, GenomePos& pos) {
    std::smatch m;
    if (std::regex_match(seq_id, m, GENOME_ID_PATTERN)) {
        pos.chr = m[1];
        pos.start = std::stoi(m[2]);
        pos.end = std::stoi(m[3]);
        return true;
    }
    return false;
}

struct ContigScore {
    int hits = 0;            // raw matched positions
    int distinct = 0;        // distinct matched mers
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
    auto [graph, kmer_counts, e_counts] =
        de_bruijn_graph(reads, k, min_count, gfa_out != nullptr, snp_out != nullptr);
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

void trim_contig(std::string& seq, const MerQueryManager& mqm, int ele_id) {
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

    int margin = 6;
    int start = std::max(0, first_hit - margin);
    int end = std::min(n, last_hit + margin);
    auto trimmed = seq.substr(start, end + mer_size - 1 - start);
    if (trimmed.size() >= 20) {
        seq = std::move(trimmed);
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
    const std::vector<std::string>& reads, int k)
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
        const double w = 1.0 / static_cast<double>(lst.size());
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
        if (kv.second < kClusterMinScore) continue;
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

    std::string best_did;
    std::string best_seq;
    double best_score = -1.0;
    size_t best_pos = 0;

    auto consider = [&](std::string contig, int cluster_start,
                        const EleRead& first_r, const EleRead& last_r, size_t pos) {
        if (contig.empty()) return;
        if (args.trim) trim_contig(contig, mqm, cur_ele);
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
        if (best_did.empty() || sc.hits > best_score ||
            (sc.hits == best_score && pos < best_pos)) {
            best_did = std::move(did);
            best_seq = std::move(contig);
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
            consider(er.seq, er.read.start, er, er, cr.pos);
            return;
        }
        std::vector<std::string> reads_for_graph = sample_cluster(ele_reads, cr.begin, cr.end);
        std::string contig = debruijn_assemble(reads_for_graph, args.kmer, args.min_count,
                                               mqm, cur_ele, nullptr, nullptr, debug_ele);
        consider(std::move(contig), ele_reads[cr.begin].read.start,
                 ele_reads[cr.begin], ele_reads[cr.end - 1], cr.pos);
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

    if (best_did.empty() && rescue_handoff_enabled()) {
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

    if (!best_did.empty()) {
        out.push_back(GenomeResult{cur_ele, std::move(best_did), std::move(best_seq)});
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
                      << "  -k, --kmer <n>        K-mer length (default: 19)\n"
                      << "  --min-count <n>       Min (k+1)-mer occurrence in graph (default: 1)\n"
                      << "  --max-reads <n>       Max reads per element/locus (default: 200)\n"
                      << "  --max-loci-gap <n>    Max gap for locus clustering (default: 5000)\n"
                      << "  --max-assign <n>      Load-time cap on (element, locus) assignments per\n"
                      << "                        k-mer, must match validate (default: 3; 0 = keep all)\n"
                      << "  -t, --threads <n>     Worker threads for genome-mode assembly (default: 4)\n"
                      << "  --trim                Trim contigs to confident k-mer region (default: on)\n"
                      << "  --no-trim             Disable trimming\n"
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

        if (!is_genome) {
            size_t idx = 0;
            while (idx < all_reads.size()) {
                int cur_ele = all_reads[idx].ele_id;

                // Reads mode: one group per element
                std::vector<std::string> reads;
                reads.reserve(args.max_reads);
                std::vector<double> rescues;   // parallel to reads: T3 rescue score
                rescues.reserve(args.max_reads);
                int cnt = 0;
                while (idx < all_reads.size() && all_reads[idx].ele_id == cur_ele && cnt < args.max_reads) {
                    const auto& r = all_reads[idx];
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
                // Skip remaining reads of this element
                while (idx < all_reads.size() && all_reads[idx].ele_id == cur_ele) ++idx;

                // Best-cluster assembly (same semantics as genome mode): split
                // the element's reads into locus-pure clusters, assemble each,
                // trim, apply the repeat-chimera guard and keep the contig with
                // the highest element-mer hit count.
                std::string contig;
                ContigScore contig_sc{0, 0};   // score of `contig` (cluster path)
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
                    std::vector<int> labels =
                        cluster_reads_by_shared_kmers(reads, args.kmer);
                    std::unordered_map<int, std::vector<int>> groups;
                    for (int i = 0; i < static_cast<int>(labels.size()); ++i)
                        groups[labels[i]].push_back(i);
                    double best_score = -1.0;
                    int second_hits = -1;
                    for (const auto& kv : groups) {
                        if (kv.second.size() < 2) continue;
                        ++dbg_ng;
                        std::vector<std::string> sub;
                        sub.reserve(kv.second.size());
                        for (int i : kv.second) sub.push_back(reads[i]);
                        std::string c = debruijn_assemble(sub, args.kmer, args.min_count,
                                                          mqm, cur_ele, snp_out_ptr,
                                                          gfa_out_ptr);
                        if (c.empty()) continue;
                        ++dbg_cand;
                        if (args.trim) trim_contig(c, mqm, cur_ele);
                        if (c.size() < 20) continue;
                        ++dbg_trim20;
                        ContigScore sc = score_contig(c, mqm, cur_ele);
                        bool repeat_chimera = sc.hits >= 10 && sc.distinct > 0 &&
                            static_cast<double>(sc.hits) / sc.distinct > 2.0;
                        if (repeat_chimera) continue;
                        ++dbg_chim;
                        double score = static_cast<double>(sc.hits);
                        if (score > best_score ||
                            (score == best_score && kv.second.size() > best_support) ||
                            (score == best_score && kv.second.size() == best_support &&
                             c.size() > contig.size())) {
                            if (best_score >= 0 &&
                                static_cast<int>(best_score) > second_hits)
                                second_hits = static_cast<int>(best_score);
                            contig = std::move(c);
                            contig_sc = sc;
                            best_score = score;
                            best_support = kv.second.size();
                        } else if (static_cast<int>(score) > second_hits) {
                            second_hits = static_cast<int>(score);
                        }
                    }
                    last_second_hits = second_hits;
                    last_total_reads = static_cast<int>(reads.size());
                    last_groups = static_cast<int>(groups.size());
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
                    if (final_sc.distinct < kMinOutputDistinct) {
                        bool keep = false;
                        if (final_sc.distinct >= kMinOutputDistinctLow &&
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
                    if (best_dist >= kMinFallbackDistinct) {
                        contig = std::move(reads[best_idx]);
                        if (args.trim) trim_contig(contig, mqm, cur_ele);
                        if (contig.size() < 20) { dbg_bail("fb_short", best_dist); continue; }
                    } else {
                        // T3 rescue handoff: validate already SW-verified these
                        // reads against the element panel, so the mer-distinct
                        // floor does not apply. Emit the best-scoring one as-is
                        // (mer-based trimming would cut it to nothing).
                        if (!rescue_handoff_enabled()) { dbg_bail("no_handoff", best_dist); continue; }
                        size_t ri = reads.size();
                        double best_rs = 0.0;
                        for (size_t i = 0; i < reads.size(); ++i) {
                            if (rescues[i] > best_rs) {
                                best_rs = rescues[i];
                                ri = i;
                            }
                        }
                        if (ri == reads.size()) { dbg_bail("no_rescue", best_dist); continue; }
                        contig = std::move(reads[ri]);
                        rescued = true;
                        if (contig.size() < 20) { dbg_bail("rc_short", (long)best_rs); continue; }
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
                        if (final_sc.distinct < kMinOutputDistinct) {
                            dbg_bail("out_floor", final_sc.distinct);
                            continue;
                        }
                    }
                }
                if (asm_debug_enabled()) {
                    std::fprintf(stderr, "ASMDBG\t%d\t%zu\t%d\t%d\t%zu\t%d\t%d\t%d\n",
                                 cur_ele, contig.size(), final_sc.hits,
                                 final_sc.distinct, best_support,
                                 last_total_reads, last_groups, last_second_hits);
                }
                results.emplace_back(std::to_string(cur_ele), contig);
            }
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
                    process_genome_element(all_reads, b.first, b.second, args, mqm, grs);
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
                                                   args, mqm, per_thread[t]);
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
        return 0;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
