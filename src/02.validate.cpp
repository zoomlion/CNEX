// MIT License
//
// Copyright (c) 2024 JiangminZheng
//
// 02.validate - Standalone C++ k-mer validator for read/genome CNE matching
// Replaces the Python hybrid prototype with a fully native implementation
// Streaming pipeline: reader thread + worker thread pool

#include <hip/panel_core.hpp>
#include <hip/validator_core.hpp>

#include <parasail.h>

#include <algorithm>
#include <atomic>
#include <charconv>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <mutex>
#include <queue>
#include <string>
#include <string_view>
#include <thread>
#include <vector>

namespace fs = std::filesystem;

// Fixed gate rows: T1 is the legacy row, T2 the built-in extra row. Both are
// always active and OR-combined; --tiers rows are appended on top.
constexpr GateTier kT1{3.0, 7, 25};
constexpr GateTier kT2{2.0, 10, 30};

struct Args {
    enum class InputType { AUTO, GENOME, FASTQ };

    std::vector<std::string> reads_files;
    std::string mers_file;
    std::string output_dir = "out";
    int64_t depth = 0;
    int threads = 4;
    int chunk_size = 10000;
    int window_size = 150;
    int step_size = 25;
    std::string pigz_path = "pigz/pigz";
    int max_diff = 2;
    double vote_frac = 0.1;
    // Load-time cap on (element, locus) assignments per k-mer; 0 = keep all.
    int max_assign = 3;
    // Extra gate rows from --tiers (appended on top of T1 and T2).
    std::vector<GateTier> extra_tiers;
    // In-validate T3 rescue: SW-verify near-miss windows/reads.
    std::string rescue_panels;        // MSA with per-element panels (leave-one-out)
    double rescue_tau = -1.0;         // normalized SW score threshold
    GateTier t3{1.5, 5, 25};          // near-miss gate row
    bool help_requested = false;
    InputType input_type = InputType::AUTO;
};

struct Read {
    std::string seq_id;
    std::string seq;
    std::string qua;
    bool is_fasta = false;
    int offset = 0;
};

// Thread-safe chunk queue with backpressure (max_size chunks in memory)

class ChunkQueue {
public:
    explicit ChunkQueue(size_t max_size = 100) : max_size_(max_size) {}

    void push(std::vector<Read> chunk) {
        std::unique_lock<std::mutex> lock(mutex_);
        cv_.wait(lock, [this] { return queue_.size() < max_size_; });
        queue_.push(std::move(chunk));
        cv_.notify_one();
    }

    bool pop(std::vector<Read>& chunk) {
        std::unique_lock<std::mutex> lock(mutex_);
        cv_.wait(lock, [this] { return !queue_.empty() || done_; });
        if (!queue_.empty()) {
            chunk = std::move(queue_.front());
            queue_.pop();
            cv_.notify_one();
            return true;
        }
        return false;
    }

    void finish() {
        std::lock_guard<std::mutex> lock(mutex_);
        done_ = true;
        cv_.notify_all();
    }

private:
    std::queue<std::vector<Read>> queue_;
    std::mutex mutex_;
    std::condition_variable cv_;
    bool done_ = false;
    size_t max_size_;
};

// Returns true if path ends with .gz

static inline bool is_gz_file(const std::string& path) {
    return path.size() >= 3 && path.compare(path.size() - 3, 3, ".gz") == 0;
}

// Buffered line reader using FILE* (from pigz pipe or fopen)

class BufferedLineReader {
public:
    explicit BufferedLineReader(const std::string& path, const std::string& pigz_path) {
        if (is_gz_file(path)) {
            std::string cmd = pigz_path + " -dc \"" + path + "\" 2>/dev/null";
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

enum class FileType { FASTQ, FASTA, UNKNOWN };

FileType detect_file_type(const std::string& path, const std::string& pigz_path) {
    FILE* fp = nullptr;
    bool is_pipe = false;

    if (is_gz_file(path)) {
        std::string cmd = pigz_path + " -dc \"" + path + "\" 2>/dev/null";
        fp = popen(cmd.c_str(), "r");
        is_pipe = true;
    } else {
        fp = fopen(path.c_str(), "r");
        is_pipe = false;
    }

    if (!fp) return FileType::UNKNOWN;
    int first = fgetc(fp);
    if (is_pipe) pclose(fp); else fclose(fp);

    if (first == '@') return FileType::FASTQ;
    if (first == '>') return FileType::FASTA;
    return FileType::UNKNOWN;
}

// ─── Reader thread: streams reads into the chunk queue ───

void reader_thread(const std::vector<std::string>& files,
                   ChunkQueue& queue,
                   int chunk_size,
                   int window_size,
                   int step_size,
                   int64_t depth,
                   const std::string& pigz_path,
                   std::atomic<int64_t>& total_reads,
                   Args::InputType input_type) {
    int64_t count = 0;

    auto flush_chunk = [&](std::vector<Read>& chunk) {
        if (!chunk.empty()) {
            queue.push(std::move(chunk));
            chunk.clear();
            chunk.reserve(chunk_size);
        }
    };

    for (const auto& file : files) {
        FileType type;
        if (input_type == Args::InputType::GENOME) {
            type = FileType::FASTA;
        } else if (input_type == Args::InputType::FASTQ) {
            type = FileType::FASTQ;
        } else {
            type = detect_file_type(file, pigz_path);
        }

        if (type == FileType::FASTQ) {
            BufferedLineReader reader(file, pigz_path);
            std::string line;
            std::vector<Read> chunk;
            chunk.reserve(chunk_size);

            while (count < depth) {
                Read read;
                if (!reader.readline(read.seq_id)) break;
                if (read.seq_id.empty()) continue;
                if (!reader.readline(read.seq)) break;
                if (!reader.readline(read.qua)) break;
                if (!reader.readline(read.qua)) break;
                read.qua.clear();
                chunk.push_back(std::move(read));
                count++;
                if (static_cast<int>(chunk.size()) >= chunk_size) {
                    flush_chunk(chunk);
                }
            }
            flush_chunk(chunk);

        } else if (type == FileType::FASTA) {
            BufferedLineReader reader(file, pigz_path);
            std::string seq_id;
            std::string seq;
            std::string line;
            std::vector<Read> chunk;
            chunk.reserve(chunk_size);
            size_t chunk_bases = 0;
            const size_t PIECE = 1u << 20;   // split big contigs into ~1MB overlapping pieces

            auto push_piece = [&](const std::string& piece, int offset, size_t total) {
                Read read;
                read.seq_id = seq_id;
                read.seq = piece;
                read.is_fasta = true;
                read.offset = offset;
                (void)total;
                chunk_bases += piece.size();
                chunk.push_back(std::move(read));
                count++;
                if (static_cast<int>(chunk.size()) >= chunk_size || chunk_bases >= (4u << 20)) {
                    flush_chunk(chunk);
                    chunk_bases = 0;
                }
            };

            auto flush_sequence = [&]() {
                if (seq_id.empty() || seq.empty() || count >= depth) return;
                const size_t L = seq.size();
                if (L <= PIECE) {
                    push_piece(seq, 0, L);
                } else {
                    const size_t overlap = static_cast<size_t>(window_size - 1);
                    for (size_t off = 0; off < L; ) {
                        size_t plen = std::min(PIECE, L - off);
                        push_piece(seq.substr(off, plen), static_cast<int>(off), L);
                        if (off + plen >= L) break;
                        off += plen - overlap;
                    }
                }
                seq.clear();
            };

            while (reader.readline(line) && count < depth) {
                if (line.empty()) continue;
                if (line[0] == '>') {
                    flush_sequence();
                    seq_id = line.substr(1);
                    auto pos = seq_id.find(' ');
                    if (pos != std::string::npos) seq_id = seq_id.substr(0, pos);
                    seq.clear();
                } else {
                    seq += line;
                }
            }
            flush_sequence();
            flush_chunk(chunk);

        } else {
            std::cerr << "Warning: unknown file type for " << file << ", skipping\n";
        }
    }

    total_reads.store(count);
    queue.finish();
}

// ─── T3 rescue verification (SIMD SW vs element panels) ───

struct VerifyContext {
    const std::unordered_map<int, std::vector<std::string>>* panels = nullptr;
    double tau = 0.9;
    parasail_matrix_t* mat = nullptr;
    std::atomic<long long>* tried = nullptr;
    std::atomic<long long>* passed = nullptr;
    std::atomic<long long>* pairs = nullptr;
};

static void clean_upper(std::string& s) {
    for (char& c : s) {
        c = (c >= 'a' && c <= 'z') ? static_cast<char>(c - 'a' + 'A') : c;
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T') c = 'A';
    }
}

// Returns the best normalized SW score over the element's panel sequences when
// it reaches tau (0.0 otherwise). The score is propagated to the output as the
// rescue tag of the read/window.
static double verify_window(const VerifyContext& vc, int ele, const std::string& w, int strand) {
    if (!vc.panels || !vc.mat) return 0.0;
    auto it = vc.panels->find(ele);
    if (it == vc.panels->end() || it->second.empty()) return 0.0;
    vc.tried->fetch_add(1, std::memory_order_relaxed);

    std::string s = w;
    clean_upper(s);
    if (strand == -1) {
        thread_local std::string _rc;
        reverse_complement(s, _rc);
        s = _rc;
    }
    double best = 0.0;
    long long np = 0;
    for (const std::string& q : it->second) {
        parasail_result_t* r = parasail_sw_striped_16(
            q.c_str(), static_cast<int>(q.size()),
            s.c_str(), static_cast<int>(s.size()), 2, 1, vc.mat);
        ++np;
        int m = std::min<int>(q.size(), s.size());
        double v = m > 0 ? static_cast<double>(r->score) / m : 0.0;
        parasail_result_free(r);
        if (v > best) best = v;
    }
    vc.pairs->fetch_add(np, std::memory_order_relaxed);
    if (best >= vc.tau) {
        vc.passed->fetch_add(1, std::memory_order_relaxed);
        return best;
    }
    return 0.0;
}

// ─── Worker thread ───

void worker_thread(int thread_id, ChunkQueue& queue, const MerQueryManager& mqm,
                   int mer_size, const ValidatorConfig& cfg,
                   const VerifyContext& vc,
                   const std::string& output_dir,
                   std::atomic<int64_t>& processed) {
    const int window_size = cfg.window_size;
    std::string output_file = output_dir + "/Assemble." + std::to_string(thread_id) + ".reads";
    FILE* fp = fopen(output_file.c_str(), "a");
    if (!fp) {
        std::cerr << "Error: cannot open output file " << output_file << "\n";
        return;
    }

    // 1MB output buffer for zero-allocation writes
    const size_t OUT_BUF_SIZE = 1 << 20;
    std::vector<char> out_buf(OUT_BUF_SIZE);
    size_t out_pos = 0;

    auto verify_fn = [&](int ele, const std::string& w, int strand) {
        return verify_window(vc, ele, w, strand);
    };

    auto flush_buffer = [&]() {
        if (out_pos > 0) {
            fwrite(out_buf.data(), 1, out_pos, fp);
            out_pos = 0;
        }
    };

    // Write one validated entry: `id strand ele seq` plus the rescue score as
    // an optional 5th column when the entry came through T3 verification.
    auto write_entry = [&](const char* sid, int strand, int ele, const char* seq,
                           double rescue_score, size_t margin) {
        if (out_pos + margin >= OUT_BUF_SIZE) flush_buffer();
        int n;
        if (rescue_score > 0.0) {
            n = snprintf(out_buf.data() + out_pos, OUT_BUF_SIZE - out_pos,
                         "%s\t%d\t%d\t%s\t%.4f\n", sid, strand, ele, seq, rescue_score);
        } else {
            n = snprintf(out_buf.data() + out_pos, OUT_BUF_SIZE - out_pos,
                         "%s\t%d\t%d\t%s\n", sid, strand, ele, seq);
        }
        if (n > 0) out_pos += static_cast<size_t>(n);
        if (out_pos > OUT_BUF_SIZE * 0.8) flush_buffer();
    };

    std::vector<Read> chunk;
    while (queue.pop(chunk)) {
        for (const auto& read : chunk) {
            if (read.is_fasta) {
                const std::string& contig = read.seq_id;
                const int base_offset = read.offset;
                const size_t margin = static_cast<size_t>(window_size) + 256;
                validate_sequence_segmented(
                    read.seq, mqm, base_offset, cfg,
                    [&](int start, int strand, int confi_id, double rescue_score) {
                        std::string sid;
                        sid.reserve(contig.size() + 24);
                        sid += contig;
                        sid += ':';
                        char num[16];
                        auto rc = std::to_chars(num, num + 15, base_offset + start + 1);
                        sid.append(num, rc.ptr - num);
                        sid += '-';
                        rc = std::to_chars(num, num + 15, base_offset + start + window_size);
                        sid.append(num, rc.ptr - num);
                        std::string wseq(read.seq, start, window_size);
                        for (char& c : wseq) {
                            if (c >= 'a' && c <= 'z') c = c - 'a' + 'A';
                        }
                        write_entry(sid.c_str(), strand, confi_id, wseq.c_str(),
                                    rescue_score, margin);
                    },
                    verify_fn);
            } else {
                double rescue_score = 0.0;
                auto [confi_id, strand] = validate_read_verified(read.seq, mqm, mer_size, -1,
                                                                  cfg, verify_fn, &rescue_score);
                if (confi_id > -1) {
                    write_entry(read.seq_id.c_str(), strand, confi_id,
                                read.seq.c_str(), rescue_score, 512);
                }
            }
        }
        processed.fetch_add(chunk.size());
    }

    flush_buffer();
    fclose(fp);
}

// ─── Progress monitor ───

void progress_monitor(std::atomic<int64_t>& processed, int64_t total) {
    auto last = std::chrono::steady_clock::now();
    int64_t last_count = 0;
    while (processed.load() < total) {
        std::this_thread::sleep_for(std::chrono::seconds(1));
        auto now = std::chrono::steady_clock::now();
        int64_t current = processed.load();
        double elapsed = std::chrono::duration<double>(now - last).count();
        double rate = (current - last_count) / elapsed;
        double pct = 100.0 * current / total;
        std::cerr << "\rProcessing reads: " << current << "/" << total
                  << " (" << pct << "%, " << rate << " reads/s)";
        last = now;
        last_count = current;
    }
    std::cerr << "\nDone.\n";
}

// ─── Command-line parsing ───

static GateTier parse_gate_row(const std::string& spec) {
    size_t p1 = spec.find(',');
    if (p1 == std::string::npos)
        throw std::runtime_error("invalid spec '" + spec + "' (expect rv,mc[,ms])");
    size_t p2 = spec.find(',', p1 + 1);
    GateTier tier;
    try {
        tier.vote_ratio = std::stod(spec.substr(0, p1));
        tier.min_c = std::stoi(spec.substr(p1 + 1,
            (p2 == std::string::npos) ? std::string::npos : p2 - p1 - 1));
        tier.min_span = (p2 == std::string::npos) ? 25 : std::stoi(spec.substr(p2 + 1));
    } catch (const std::exception&) {
        throw std::runtime_error("invalid spec '" + spec + "' (expect rv,mc[,ms])");
    }
    if (tier.vote_ratio <= 0 || tier.min_c <= 0 || tier.min_span <= 0)
        throw std::runtime_error("invalid spec '" + spec + "' (expect rv,mc[,ms])");
    return tier;
}

// Parse a space-separated list of rows, each "rv,mc[,ms]".
static std::vector<GateTier> parse_tiers(const std::string& spec) {
    std::vector<GateTier> out;
    size_t i = 0;
    while (i < spec.size()) {
        while (i < spec.size() && (spec[i] == ' ' || spec[i] == '\t')) ++i;
        size_t j = i;
        while (j < spec.size() && spec[j] != ' ' && spec[j] != '\t') ++j;
        if (j > i) out.push_back(parse_gate_row(spec.substr(i, j - i)));
        i = j;
    }
    if (out.empty())
        throw std::runtime_error("--tiers: empty tier list");
    return out;
}

Args parse_args(int argc, char* argv[]) {
    Args args;
    int i = 1;
    while (i < argc) {
        std::string arg(argv[i]);
        if (arg == "-h" || arg == "--help") {
            args.help_requested = true;
            return args;
        } else if (arg == "--mers" || arg == "-m") {
            if (++i >= argc) throw std::runtime_error("Missing value for --mers");
            args.mers_file = argv[i];
        } else if (arg == "--depth" || arg == "-d") {
            if (++i >= argc) throw std::runtime_error("Missing value for --depth");
            args.depth = std::stoll(argv[i]);
        } else if (arg == "-t" || arg == "--threads") {
            if (++i >= argc) throw std::runtime_error("Missing value for --threads");
            args.threads = std::stoi(argv[i]);
        } else if (arg == "--output_dir" || arg == "-o") {
            if (++i >= argc) throw std::runtime_error("Missing value for --output_dir");
            args.output_dir = argv[i];
        } else if (arg == "--chunk_size" || arg == "-c") {
            if (++i >= argc) throw std::runtime_error("Missing value for --chunk_size");
            args.chunk_size = std::stoi(argv[i]);
        } else if (arg == "--window_size" || arg == "-w") {
            if (++i >= argc) throw std::runtime_error("Missing value for --window_size");
            args.window_size = std::stoi(argv[i]);
        } else if (arg == "--step_size" || arg == "-s") {
            if (++i >= argc) throw std::runtime_error("Missing value for --step_size");
            args.step_size = std::stoi(argv[i]);
        } else if (arg == "--type") {
            if (++i >= argc) throw std::runtime_error("Missing value for --type");
            std::string type = argv[i];
            if (type == "genome") args.input_type = Args::InputType::GENOME;
            else if (type == "fastq") args.input_type = Args::InputType::FASTQ;
            else throw std::runtime_error("Unknown --type: " + type + " (use genome|fastq)");
        } else if (arg == "--pigz") {
            if (++i >= argc) throw std::runtime_error("Missing value for --pigz");
            args.pigz_path = argv[i];
        } else if (arg == "--max-diff") {
            if (++i >= argc) throw std::runtime_error("Missing value for --max-diff");
            args.max_diff = std::stoi(argv[i]);
        } else if (arg == "--vote-frac") {
            if (++i >= argc) throw std::runtime_error("Missing value for --vote-frac");
            args.vote_frac = std::stod(argv[i]);
        } else if (arg == "--max-assign") {
            if (++i >= argc) throw std::runtime_error("Missing value for --max-assign");
            args.max_assign = std::stoi(argv[i]);
        } else if (arg == "--rescue-panels") {
            if (++i >= argc) throw std::runtime_error("Missing value for --rescue-panels");
            args.rescue_panels = argv[i];
        } else if (arg == "--rescue-tau") {
            if (++i >= argc) throw std::runtime_error("Missing value for --rescue-tau");
            args.rescue_tau = std::stod(argv[i]);
        } else if (arg == "--t3") {
            if (++i >= argc) throw std::runtime_error("Missing value for --t3");
            args.t3 = parse_gate_row(argv[i]);
        } else if (arg == "--tiers") {
            if (++i >= argc) throw std::runtime_error("Missing value for --tiers");
            for (auto& t : parse_tiers(argv[i])) args.extra_tiers.push_back(t);
        } else if (arg[0] != '-') {
            // Heuristic: an unquoted --tiers list arrives as extra tokens such
            // as "1.5,12,40" and would silently be treated as input files.
            if (arg.find(',') != std::string::npos &&
                arg.find_first_not_of("0123456789.,") == std::string::npos) {
                throw std::runtime_error(
                    "unexpected argument '" + arg +
                    "' - quote the whole --tiers list (e.g. --tiers \"2.0,10,30 1.5,12,40\")");
            }
            args.reads_files.push_back(arg);
        } else {
            throw std::runtime_error("Unknown option: " + arg);
        }
        ++i;
    }
    return args;
}

// ─── Gate matrix ───

constexpr size_t kMaxGateTiers = 6;

// T1 and T2 are always active; --tiers rows are appended on top of them.
static ValidatorConfig build_validator_config(const Args& args) {
    ValidatorConfig cfg;
    cfg.window_size = args.window_size;
    cfg.step_size = args.step_size;
    cfg.vote_frac = args.vote_frac;
    cfg.max_diff = args.max_diff;
    cfg.tiers.clear();
    cfg.tiers.push_back(kT1);
    cfg.tiers.push_back(kT2);
    for (const auto& tier : args.extra_tiers) {
        cfg.tiers.push_back(tier);
    }
    if (cfg.tiers.size() > kMaxGateTiers) {
        throw std::runtime_error("too many gate rows: T1 + T2 + --tiers must be <= 6");
    }
    // Sort extra rows by vote_ratio descending, then enforce monotonicity: a
    // lower vote_ratio must not lower the required chain length.
    std::sort(cfg.tiers.begin() + 1, cfg.tiers.end(),
              [](const GateTier& a, const GateTier& b) {
                  return a.vote_ratio > b.vote_ratio;
              });
    for (size_t t = 1; t < cfg.tiers.size(); ++t) {
        if (cfg.tiers[t].vote_ratio < cfg.tiers[t - 1].vote_ratio &&
            cfg.tiers[t].min_c < cfg.tiers[t - 1].min_c) {
            throw std::runtime_error(
                "gate rows must be monotone: lower vote_ratio requires >= min_c");
        }
    }
    return cfg;
}

// ─── Main ───

int main(int argc, char* argv[]) {
    try {
        Args args = parse_args(argc, argv);

        if (args.help_requested || args.reads_files.empty()) {
            std::cerr << "Usage: " << argv[0]
                      << " <input_files...> --mers <mers_table> [options]\n"
                      << "\nOptions:\n"
                      << "  --mers <file>         TSV mers table (required)\n"
                      << "  --type genome|fastq   Input type: genome (sliding-window FASTA) or\n"
                      << "                        fastq (reads). Auto-detected if omitted.\n"
                      << "  --depth <n>           Max reads/pseudo-reads to process (default: 0=unlimited)\n"
                      << "  -t, --threads <n>     Number of threads (default: 4)\n"
                      << "  --output_dir <dir>    Output directory (default: out)\n"
                      << "  --chunk_size <n>      Reads per chunk (default: 10000)\n"
                      << "  --window_size <n>     Genome sliding window size (default: 150)\n"
                      << "  --step_size <n>       Genome sliding window step (default: 25, ~6X coverage)\n"
                      << "  --pigz <path>         pigz binary path (default: pigz/pigz)\n"
                      << "  --max-diff <n>        chain gap tolerance, all rows (default: 2)\n"
                      << "  --vote-frac <f>       min fraction of k-mers for top CNE vote, all rows\n"
                      << "                        (default: 0.1)\n"
                      << "  --max-assign <n>      Load-time cap: k-mers with more than n (element, locus)\n"
                      << "                        assignments are dropped (default: 3; 0 = keep all)\n"
                      << "  --tiers \"<rv,mc[,ms]> <rv,mc[,ms]> ...\"\n"
                      << "                        space-separated extra OR gate rows appended on top\n"
                      << "                        of the fixed T1 (3.0,7,25) and T2 (2.0,10,30) rows;\n"
                      << "                        quote the whole list. ms defaults to 25. Rows are\n"
                      << "                        sorted by vote_ratio (desc); a row with lower\n"
                      << "                        vote_ratio must not have a lower min_c than the row\n"
                      << "                        above. Max 6 rows in total.\n"
                      << "  --rescue-panels <msa> In-validate T3 rescue: per-element panels from the\n"
                      << "                        leave-one-out MSA (top-3 longest non-ref sequences)\n"
                      << "  --rescue-tau <f>      Normalized SW score threshold for the T3 rescue\n"
                      << "  --t3 \"<rv,mc[,ms]>\"   Near-miss gate row (default: 1.5,5,25)\n";
            return 1;
        }

        if (args.mers_file.empty()) {
            std::cerr << "Error: --mers is required\n";
            return 1;
        }

        // No depth limit: process all reads/windows
        if (args.depth <= 0) {
            args.depth = std::numeric_limits<int64_t>::max();
        }

        // Validate CLI-derived gate matrix before touching disk or loading data
        ValidatorConfig cfg = build_validator_config(args);

        // In-validate T3 rescue: load per-element panels and prepare the SIMD
        // SW verification context (shared read-only across workers).
        std::unordered_map<int, std::vector<std::string>> panels;
        VerifyContext vc;
        std::atomic<long long> v_tried{0}, v_passed{0}, v_pairs{0};
        if (!args.rescue_panels.empty()) {
            if (args.rescue_tau < 0.0) {
                std::cerr << "Error: --rescue-panels requires --rescue-tau\n";
                return 1;
            }
            std::cerr << "Loading rescue panels from " << args.rescue_panels << " ...\n";
            panel_core::load_panels(args.rescue_panels, 3, panels);
            size_t nseq = 0;
            for (const auto& kv : panels) nseq += kv.second.size();
            std::cerr << "  panels: " << panels.size() << " elements, " << nseq
                      << " sequences; T3=" << args.t3.vote_ratio << ","
                      << args.t3.min_c << "," << args.t3.min_span
                      << " tau=" << args.rescue_tau << "\n";
            cfg.verify.enabled = true;
            cfg.verify.t3 = args.t3;
            cfg.verify.tau = args.rescue_tau;
            vc.panels = &panels;
            vc.tau = args.rescue_tau;
            vc.mat = parasail_matrix_create("ACGT", 2, -1);
            if (!vc.mat) throw std::runtime_error("parasail_matrix_create failed");
            vc.tried = &v_tried;
            vc.passed = &v_passed;
            vc.pairs = &v_pairs;
        }

        // Load mers table
        std::cerr << "Loading mers table from " << args.mers_file << " ...\n";
        MerQueryManager mqm;
        mqm.load_from_file(args.mers_file, args.max_assign);
        int mer_size = mqm.get_mer_size();
        std::cerr << "  Loaded " << mqm.size() << " mers, k=" << mer_size << "\n";

        // Build the direct-index accelerator for genome-mode window scans
        // (k <= 13). Reads mode does not use it, so skip it there.
        bool genome_input = (args.input_type == Args::InputType::GENOME);
        if (!genome_input && !args.reads_files.empty() &&
            args.input_type == Args::InputType::AUTO &&
            detect_file_type(args.reads_files[0], args.pigz_path) == FileType::FASTA) {
            genome_input = true;
        }
        if (genome_input && mer_size >= 7) {
            mqm.build_direct_index();
            if (!mqm.has_direct_index())
                std::cerr << "  (direct index unavailable for k=" << mer_size
                          << ", using hash lookups)\n";
        }

        // Create output directory
        if (fs::exists(args.output_dir)) {
            fs::remove_all(args.output_dir);
        }
        fs::create_directories(args.output_dir);

        // Initialize output files
        for (int i = 0; i < args.threads; ++i) {
            std::string path = args.output_dir + "/Assemble." + std::to_string(i) + ".reads";
            std::ofstream ofs(path, std::ios::trunc);
            if (!ofs) {
                throw std::runtime_error("Cannot create output file: " + path);
            }
        }

        // Pipeline: workers start first, then reader streams in chunks
        ChunkQueue queue(16);
        std::atomic<int64_t> processed{0};
        std::atomic<int64_t> total_reads{0};
        std::vector<std::thread> workers;

        for (int i = 0; i < args.threads; ++i) {
            workers.emplace_back(worker_thread, i, std::ref(queue),
                                 std::cref(mqm), mer_size,
                                 std::cref(cfg),
                                 std::cref(vc),
                                 std::cref(args.output_dir),
                                 std::ref(processed));
        }

        // Start reader thread (reads and pushes chunks)
        std::cerr << "Reading and processing reads ...\n";
        std::thread reader(reader_thread, std::cref(args.reads_files), std::ref(queue),
                           args.chunk_size, args.window_size, args.step_size,
                           args.depth, std::cref(args.pigz_path), std::ref(total_reads),
                           args.input_type);

        // Wait for reader to finish (all chunks in queue)
        reader.join();
        int64_t n_total = total_reads.load();
        std::cerr << "  Total reads to process: " << n_total << "\n";

        // Now start progress monitor (we know total)
        std::thread monitor(progress_monitor, std::ref(processed), n_total);

        // Wait for workers to finish processing queue
        for (auto& t : workers) {
            t.join();
        }

        // Wait for monitor
        processed.store(n_total);
        monitor.join();

        if (vc.mat) {
            parasail_matrix_free(vc.mat);
            std::cerr << "  T3 rescue: verified=" << v_tried
                      << ", passed=" << v_passed << ", sw_pairs=" << v_pairs << "\n";
        }

        if (validate_profile_enabled()) {
            std::cerr << "[profile] windows=" << g_val_profile.windows
                      << " gathers=" << g_val_profile.gathers
                      << " emits=" << g_val_profile.emits
                      << "\n[profile] pushed=" << g_val_profile.pushed
                      << " tally_calls=" << g_val_profile.tally_calls
                      << " chain_calls=" << g_val_profile.chain_calls
                      << " prep_ns=" << g_val_profile.prep_ns
                      << " win_ns=" << g_val_profile.win_ns
                      << "\n[profile] gather_ns=" << g_val_profile.gather_ns
                      << " tally_ns=" << g_val_profile.tally_ns << "\n";
        }

        std::cerr << "Done. Output written to " << args.output_dir << "\n";
        return 0;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
