#include <hip/debruijn_core.hpp>
#include <hip/validator_core.hpp>

#include <tsl/robin_set.h>
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <limits>
#include <map>
#include <regex>
#include <string>
#include <vector>

namespace fs = std::filesystem;

struct Args {
    std::string input_dir;
    std::string mers_file;
    std::string output = "assembled.fasta";
    int kmer = 21;
    int min_c = 3;
    int min_count = 2;
    int max_reads = 200;
    int max_loci_gap = 5000;
    bool trim = true;
    bool trim_set = false;
    bool no_trim_set = false;
    double repeat_ratio = 1.2;
    bool output_snp = false;
    bool help_requested = false;
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

int score_contig(const std::string& seq, const MerQueryManager& mqm, int ele_id) {
    int sc = 0;
    int k = mqm.get_mer_size();
    int mer_size = mqm.get_mer_size();
    for (size_t i = 0; i + mer_size <= seq.size(); ++i) {
        try {
            auto code = dna_encoder(seq.substr(i, mer_size));
            auto it = mqm.compressed_mer_query.find(code);
            if (it != mqm.compressed_mer_query.end() && it->second.first == ele_id) {
                ++sc;
            }
        } catch (...) {
            continue;
        }
    }
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
                if (it != mqm.compressed_mer_query.end() && it->second.first == ele_id) {
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
                               std::ofstream* gfa_out = nullptr)
{
    auto [graph, kmer_counts, e_counts] = de_bruijn_graph(reads, k, min_count);
    if (graph.empty()) return "";

    std::unordered_map<std::string, int> node_scores;
    auto ns = build_node_scores(graph, mqm, ele_id, k);
    if (!ns.empty()) node_scores = std::move(ns);

    std::vector<std::string> contig_path;
    std::string contig = assemble_sequence(graph, node_scores, &contig_path);

    bool has_var = false;
    if (snp_out && !contig_path.empty())
        has_var = scan_snps(graph, e_counts, contig_path, contig, ele_id, k, *snp_out);

    if (has_var && gfa_out && !contig_path.empty())
        export_path_gfa(graph, node_scores, contig_path, ele_id, k, *gfa_out);

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
            if (it != mqm.compressed_mer_query.end() && it->second.first == ele_id) {
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
    auto [confi_id, _] = validate_read(seq, mqm, mer_size, min_c);
    return confi_id == ele_id;
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
                      << "  -k, --kmer <n>        K-mer length (default: 21)\n"
                      << "  --min-c <n>           Min k-mer matches for validation (default: 3)\n"
                      << "  --min-count <n>       Min (k+1)-mer occurrence in graph (default: 2)\n"
                      << "  --max-reads <n>       Max reads per element/locus (default: 200)\n"
                      << "  --max-loci-gap <n>    Max gap for locus clustering (default: 5000)\n"
                      << "  --trim                Trim contigs to confident k-mer region (default: on)\n"
                      << "  --no-trim             Disable trimming\n"
                      << "  --repeat-ratio <f>    Max total/unique mers per read (default: 1.2)\n"
                      << "  --snp                 Scan for candidate SNPs from validated reads\n";
            return 1;
        }

        std::cerr << "Loading mers table from " << args.mers_file << " ...\n";
        MerQueryManager mqm;
        mqm.load_from_file(args.mers_file);
        int mer_size = mqm.get_mer_size();
        std::cerr << "  Loaded " << mqm.size() << " mers, k=" << mer_size << "\n";

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
        struct BestEntry { std::string display_id; std::string seq; int score; };
        std::map<int, BestEntry> best;

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

        size_t idx = 0;
        while (idx < all_reads.size()) {
            // Start of a new element
            int cur_ele = all_reads[idx].ele_id;

            if (!is_genome) {
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
                                if (it != mqm.compressed_mer_query.end() && it->second.first == cur_ele) {
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

                if (args.trim) trim_contig(contig, mqm, cur_ele);
                if (contig.size() < 20) continue;

                results.emplace_back(std::to_string(cur_ele), contig);
            } else {
                // Genome mode: collect all reads for this element
                std::vector<std::pair<Read, std::string>> ele_reads; // (read, local_seq)
                ele_reads.reserve(args.max_reads);
                int cnt = 0;
                while (idx < all_reads.size() && all_reads[idx].ele_id == cur_ele && cnt < args.max_reads) {
                    auto& r = all_reads[idx];
                    std::string local_seq = r.seq;
                    for (char& c : local_seq) {
                        if (c >= 'a' && c <= 'z') c = c - 'a' + 'A';
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
                                if (it != mqm.compressed_mer_query.end() && it->second.first == cur_ele) {
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
                while (idx < all_reads.size() && all_reads[idx].ele_id == cur_ele) ++idx;

                // Cluster by locus
                auto sort_key = [](const std::pair<Read, std::string>& a) {
                    return std::make_tuple(a.first.chr, a.first.strand, a.first.start);
                };
                std::sort(ele_reads.begin(), ele_reads.end(),
                    [](const auto& a, const auto& b) {
                        if (a.first.chr != b.first.chr) return a.first.chr < b.first.chr;
                        if (a.first.strand != b.first.strand) return a.first.strand > b.first.strand;
                        return a.first.start < b.first.start;
                    });

                // Per-locus assembly
                std::vector<std::pair<Read, std::string>> cluster;
                int cluster_start = ele_reads[0].first.start;
                std::string best_display_id;

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
                        // Assemble current cluster
                        std::vector<std::string> reads_for_graph;
                        reads_for_graph.reserve(cluster.size());
                        for (const auto& c : cluster) {
                            reads_for_graph.push_back(c.second);
                        }

                        std::string contig = debruijn_assemble(reads_for_graph, args.kmer, args.min_count, mqm, cur_ele, snp_out_ptr, gfa_out_ptr);
                        if (!contig.empty()) {
                            if (args.trim) trim_contig(contig, mqm, cur_ele);
                            if (contig.size() >= 20) {
                                int sc = score_contig(contig, mqm, cur_ele);
                                char s = cluster[0].first.strand == 1 ? '+' : '-';
                                std::string did = std::to_string(cur_ele) + "." +
                                    cluster[0].first.chr + ":" +
                                    std::to_string(cluster_start) + "-" +
                                    std::to_string(cluster.back().first.end) + "(" + s + ")";
                                auto it = best.find(cur_ele);
                                if (it == best.end() || sc > it->second.score) {
                                    best[cur_ele] = BestEntry{did, contig, sc};
                                }
                            }
                        }

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

                    std::string contig = debruijn_assemble(reads_for_graph, args.kmer, args.min_count, mqm, cur_ele, snp_out_ptr, gfa_out_ptr);
                    if (!contig.empty()) {
                        if (args.trim) trim_contig(contig, mqm, cur_ele);
                        if (contig.size() >= 20) {
                            int sc = score_contig(contig, mqm, cur_ele);
                            char s = cluster[0].first.strand == 1 ? '+' : '-';
                            std::string did = std::to_string(cur_ele) + "." +
                                cluster[0].first.chr + ":" +
                                std::to_string(cluster_start) + "-" +
                                std::to_string(cluster.back().first.end) + "(" + s + ")";
                            auto it = best.find(cur_ele);
                            if (it == best.end() || sc > it->second.score) {
                                best[cur_ele] = BestEntry{did, contig, sc};
                            }
                        }
                    }
                }
            }
        }

        // For genome mode: collect best results
        if (is_genome) {
            for (const auto& [eid, entry] : best) {
                results.emplace_back(entry.display_id, entry.seq);
            }
            std::sort(results.begin(), results.end(),
                [](const auto& a, const auto& b) {
                    int id_a = std::stoi(a.first.substr(0, a.first.find('.')));
                    int id_b = std::stoi(b.first.substr(0, b.first.find('.')));
                    return id_a < id_b;
                });
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
