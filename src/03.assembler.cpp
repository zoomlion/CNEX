#include <hip/debruijn_core.hpp>
#include <hip/validator_core.hpp>

#include <tsl/robin_set.h>
#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdio>
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
    int min_c = 3;
    int min_count = 2;
    int max_reads = 200;
    int max_loci_gap = 5000;
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
};

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
        } else if (arg == "--min-c") {
            if (++i >= argc) throw std::runtime_error("Missing value for --min-c");
            args.min_c = std::stoi(argv[i]);
        } else if (arg == "--min-count") {
            if (++i >= argc) throw std::runtime_error("Missing value for --min-count");
            args.min_count = std::stoi(argv[i]);
        } else if (arg == "--max-reads") {
            if (++i >= argc) throw std::runtime_error("Missing value for --max-reads");
            args.max_reads = std::stoi(argv[i]);
        } else if (arg == "--max-loci-gap") {
            if (++i >= argc) throw std::runtime_error("Missing value for --max-loci-gap");
            args.max_loci_gap = std::stoi(argv[i]);
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
        try {
            auto code = dna_encoder(seq.substr(i, mer_size));
            auto it = mqm.compressed_mer_query.find(code);
            if (it != mqm.compressed_mer_query.end() && mqm.has_id(it->second, ele_id)) {
                ++sc.hits;
                seen.insert(code);
            }
        } catch (...) {
            continue;
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
            try {
                auto code = dna_encoder(node.substr(i, mer_size));
                auto it = mqm.compressed_mer_query.find(code);
                if (it != mqm.compressed_mer_query.end() && mqm.has_id(it->second, ele_id)) {
                    ++score;
                }
            } catch (...) {
                continue;
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
    auto [graph, kmer_counts, e_counts] = de_bruijn_graph(reads, k, min_count);
    if (graph.empty()) return "";

    std::unordered_map<std::string, int> node_scores;
    auto ns = build_node_scores(graph, mqm, ele_id, k);
    if (!ns.empty()) node_scores = std::move(ns);

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
    int n = static_cast<int>(seq.size()) - 12;

    int first_hit = -1, last_hit = -1;
    for (int i = 0; i < n; ++i) {
        try {
            auto code = dna_encoder(seq.substr(i, 13));
            auto it = mqm.compressed_mer_query.find(code);
            if (it != mqm.compressed_mer_query.end() && mqm.has_id(it->second, ele_id)) {
                if (first_hit < 0) first_hit = i;
                last_hit = i;
            }
        } catch (...) {
            continue;
        }
    }

    if (first_hit < 0) return;

    int margin = 6;
    int start = std::max(0, first_hit - margin);
    int end = std::min(n, last_hit + margin);
    auto trimmed = seq.substr(start, end + 12 - start);
    if (trimmed.size() >= 20) {
        seq = std::move(trimmed);
    }
}

bool filter_contig(const std::string& seq, const MerQueryManager& mqm, int mer_size,
                   int ele_id, int min_c)
{
    if ((int)seq.size() < mer_size + min_c) return true;
    ValidatorConfig cfg;
    cfg.tiers[0].min_c = min_c;
    auto [confi_id, _] = validate_read(seq, mqm, mer_size, -1, cfg);
    return confi_id == ele_id;
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

void sw_align(const std::string& ref, const std::string& seq,
              int& r_start, int& r_end, int& score) {
    int m = (int)ref.size(), n = (int)seq.size();
    const int match_score = 5, mismatch_score = -4;
    const int gap_open = 15, gap_extend = 1;

    std::vector<std::vector<int>> M(m + 1, std::vector<int>(n + 1, 0));
    std::vector<std::vector<int>> Ix(m + 1, std::vector<int>(n + 1, 0));
    std::vector<std::vector<int>> Iy(m + 1, std::vector<int>(n + 1, 0));

    score = 0;
    int max_i = 0, max_j = 0;

    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            int d = (ref[i - 1] == seq[j - 1]) ? match_score : mismatch_score;
            M[i][j] = std::max({M[i - 1][j - 1] + d, Ix[i - 1][j - 1] + d,
                                Iy[i - 1][j - 1] + d, 0});
            Ix[i][j] = std::max({M[i - 1][j] - gap_open, Ix[i - 1][j] - gap_extend, 0});
            Iy[i][j] = std::max({M[i][j - 1] - gap_open, Iy[i][j - 1] - gap_extend, 0});
            int cur = std::max({M[i][j], Ix[i][j], Iy[i][j]});
            if (cur > score) { score = cur; max_i = i; max_j = j; }
        }
    }
    if (score == 0) { r_start = 0; r_end = n; return; }

    // Traceback from (max_i, max_j) to find start position in seq
    int i = max_i, j = max_j;
    int cur = std::max({M[i][j], Ix[i][j], Iy[i][j]});
    while (i > 0 && j > 0 && cur > 0) {
        if (cur == M[i][j]) {
            int d = (ref[i - 1] == seq[j - 1]) ? match_score : mismatch_score;
            int prev = std::max({M[i - 1][j - 1], Ix[i - 1][j - 1], Iy[i - 1][j - 1]});
            if (M[i][j] == prev + d) cur = M[i - 1][j - 1];
            else cur = std::max({M[i - 1][j - 1], Ix[i - 1][j - 1], Iy[i - 1][j - 1]});
            --i; --j;
        } else if (cur == Ix[i][j]) {
            cur = Ix[i - 1][j];
            --i;
        } else {
            cur = Iy[i][j - 1];
            --j;
        }
        if (cur <= 0) break;
    }
    r_start = j;
    r_end = max_j;
}

void trim_contig_sw(std::string& seq, const std::string& ref,
                    const MerQueryManager& mqm, int ele_id) {
    if (seq.size() < 26 || ref.empty()) return;
    int mer_size = mqm.get_mer_size();
    int n = static_cast<int>(seq.size()) - mer_size + 1;

    // Phase 1: 13-mer扫描确认核心区
    int first_hit = -1, last_hit = -1;
    for (int i = 0; i < n; ++i) {
        try {
            auto code = dna_encoder(seq.substr(i, 13));
            auto it = mqm.compressed_mer_query.find(code);
            if (it != mqm.compressed_mer_query.end() && mqm.has_id(it->second, ele_id)) {
                if (first_hit < 0) first_hit = i;
                last_hit = i;
            }
        } catch (...) { continue; }
    }
    if (first_hit < 0) return;

    // Phase 2: SW 扩展边界
    const int extend = 100;
    int L = std::max(0, first_hit - extend);
    int R = std::min((int)seq.size(), last_hit + mer_size + extend);
    if (R - L < 20) return;

    std::string candidate = seq.substr(L, R - L);
    int r_start, r_end, score;
    sw_align(ref, candidate, r_start, r_end, score);
    if (score == 0) return;

    int sw_start = L + r_start;
    int sw_end = L + r_end + 1;

    // 保底下界：只扩展，不内缩
    int old_start = first_hit - 6;
    int old_end = last_hit + mer_size + 6;
    int new_start = std::min(sw_start, old_start);
    int new_end = std::max(sw_end, old_end);

    new_start = std::max(0, new_start);
    new_end = std::min((int)seq.size(), new_end);

    if (new_end - new_start < 20) return;
    seq = seq.substr(new_start, new_end - new_start);
}


struct GenomeResult {
    int ele_id;
    std::string display_id;
    std::string seq;
};

static void process_genome_element(
    const std::vector<Read>& all_reads, size_t begin, size_t end,
    const Args& args, const MerQueryManager& mqm,
    const std::unordered_map<int, std::string>& cne_refs,
    std::vector<GenomeResult>& out)
{
    const int cur_ele = all_reads[begin].ele_id;
    std::vector<std::pair<Read, std::string>> ele_reads;
    ele_reads.reserve(args.max_reads);
    int cnt = 0;
    size_t idx = begin;
    while (idx < end && cnt < args.max_reads) {
        const auto& r = all_reads[idx];
        std::string local_seq = r.seq;
        for (char& c : local_seq) {
            if (c >= 'a' && c <= 'z') c = c - 'a' + 'A';
        }
        for (char& c : local_seq) {
            if (c != 'A' && c != 'C' && c != 'G' && c != 'T') c = 'N';
        }
        if (r.strand == -1) {
            thread_local std::string _rc;
            reverse_complement(local_seq, _rc);
            local_seq = _rc;
        }
        // Repeat ratio check for genome mode
        {
            tsl::robin_set<uint32_t> seen;
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
                ++cnt;
                continue;
            }
        }
        ele_reads.emplace_back(r, std::move(local_seq));
        ++idx;
        ++cnt;
    }
    if (ele_reads.empty()) return;

    // Cluster by locus
    std::sort(ele_reads.begin(), ele_reads.end(),
        [](const auto& a, const auto& b) {
            if (a.first.chr != b.first.chr) return a.first.chr < b.first.chr;
            if (a.first.strand != b.first.strand) return a.first.strand > b.first.strand;
            return a.first.start < b.first.start;
        });

    std::string best_did;
    std::string best_seq;
    double best_score = -1.0;

    auto consider = [&](std::string contig, int cluster_start,
                        const std::pair<Read, std::string>& first_r,
                        const std::pair<Read, std::string>& last_r) {
        if (contig.empty()) return;
        if (args.trim) {
            auto it = cne_refs.find(cur_ele);
            if (it != cne_refs.end()) {
                trim_contig_sw(contig, it->second, mqm, cur_ele);
            } else {
                trim_contig(contig, mqm, cur_ele);
            }
        }
        if (contig.size() < 20) return;
        ContigScore sc = score_contig(contig, mqm, cur_ele);
        bool repeat_chimera = sc.hits >= 10 && sc.distinct > 0 &&
            static_cast<double>(sc.hits) / sc.distinct > 2.0;
        if (repeat_chimera) return;
        char s = first_r.first.strand == 1 ? '+' : '-';
        std::string did = std::to_string(cur_ele) + "." +
            first_r.first.chr + ":" +
            std::to_string(cluster_start) + "-" +
            std::to_string(last_r.first.end) + "(" + s + ")";
        if (best_did.empty() || sc.hits > best_score) {
            best_did = std::move(did);
            best_seq = std::move(contig);
            best_score = static_cast<double>(sc.hits);
        }
    };

    // Per-locus assembly
    std::vector<std::pair<Read, std::string>> cluster;
    int cluster_start = ele_reads[0].first.start;

    for (size_t j = 0; j < ele_reads.size(); ++j) {
        const auto& r = ele_reads[j];
        bool new_cluster = false;
        if (!cluster.empty()) {
            int gap = r.first.start - cluster.back().first.end;
            if (gap > args.max_loci_gap ||
                r.first.chr != cluster.back().first.chr ||
                r.first.strand != cluster.back().first.strand)
            {
                new_cluster = true;
            }
        }

        if (new_cluster) {
            std::vector<std::string> reads_for_graph;
            reads_for_graph.reserve(cluster.size());
            for (const auto& c : cluster) {
                reads_for_graph.push_back(c.second);
            }
            std::string contig = debruijn_assemble(reads_for_graph, args.kmer, args.min_count,
                                                   mqm, cur_ele, nullptr, nullptr,
                                                   args.debug_eles.count(cur_ele) > 0);
            consider(std::move(contig), cluster_start, cluster.front(), cluster.back());
            cluster.clear();
            cluster_start = r.first.start;
        }
        cluster.push_back(r);
    }

    // Last cluster
    if (!cluster.empty()) {
        std::vector<std::string> reads_for_graph;
        reads_for_graph.reserve(cluster.size());
        for (const auto& c : cluster) {
            reads_for_graph.push_back(c.second);
        }
        std::string contig = debruijn_assemble(reads_for_graph, args.kmer, args.min_count,
                                               mqm, cur_ele, nullptr, nullptr,
                                               args.debug_eles.count(cur_ele) > 0);
        consider(std::move(contig), cluster_start, cluster.front(), cluster.back());
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
                      << "  --min-c <n>           Min k-mer matches for validation (default: 3)\n"
                      << "  --min-count <n>       Min (k+1)-mer occurrence in graph (default: 2)\n"
                      << "  --max-reads <n>       Max reads per element/locus (default: 200)\n"
                      << "  --max-loci-gap <n>    Max gap for locus clustering (default: 5000)\n"
                      << "  -t, --threads <n>     Worker threads for genome-mode assembly (default: 4)\n"
                      << "  --trim                Trim contigs to confident k-mer region (default: on)\n"
                      << "  --no-trim             Disable trimming\n"
                      << "  --repeat-ratio <f>    Max total/unique mers per read (default: 1.2)\n"
                       << "  --snp                 Scan for candidate SNPs from validated reads\n"
                       << "  --cne-ref <file>      CNE reference FASTA for SW-based trimming (ref.fa)\n";
            return 1;
        }

        std::cerr << "Loading mers table from " << args.mers_file << " ...\n";
        MerQueryManager mqm;
        mqm.load_from_file(args.mers_file);
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
                std::string seq = line.substr(t3 + 1);

                if (seq.empty()) continue;

                Read read;
                read.seq_id = seq_id;
                read.strand = strand;
                read.ele_id = ele_id;
                read.seq = seq;

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
                int cnt = 0;
                while (idx < all_reads.size() && all_reads[idx].ele_id == cur_ele && cnt < args.max_reads) {
                    const auto& r = all_reads[idx];
                    std::string local_seq = r.seq;
                    for (char& c : local_seq) {
                        if (c >= 'a' && c <= 'z') c = c - 'a' + 'A';
                    }
                    for (char& c : local_seq) {
                        if (c != 'A' && c != 'C' && c != 'G' && c != 'T') c = 'N';
                    }
                    if (r.strand == -1) {
                        thread_local std::string _rc;
                        reverse_complement(local_seq, _rc);
                        local_seq = _rc;
                    }
                    // Repeat ratio check
                    {
                        tsl::robin_set<uint32_t> seen;
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
                            ++cnt;
                            continue;
                        }
                    }
                    reads.push_back(std::move(local_seq));
                    ++idx;
                    ++cnt;
                }
                // Skip remaining reads of this element
                while (idx < all_reads.size() && all_reads[idx].ele_id == cur_ele) ++idx;

                std::string contig = debruijn_assemble(reads, args.kmer, args.min_count, mqm, cur_ele, snp_out_ptr, gfa_out_ptr);
                if (contig.empty()) continue;

                if (args.trim) {
                    auto it = cne_refs.find(cur_ele);
                    if (it != cne_refs.end()) {
                        trim_contig_sw(contig, it->second, mqm, cur_ele);
                    } else {
                        trim_contig(contig, mqm, cur_ele);
                    }
                }
                if (contig.size() < 20) continue;

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
                    process_genome_element(all_reads, b.first, b.second, args, mqm, cne_refs, grs);
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
                                                   args, mqm, cne_refs, per_thread[t]);
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
