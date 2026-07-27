// MIT License
//
// Copyright (c) 2024 JiangminZheng
//
// 02.validate - Standalone C++ k-mer validator for read/genome CNE matching
// Replaces the Python hybrid prototype with a fully native implementation
// Streaming pipeline: reader thread + worker thread pool

#include <hip/validator_core.hpp>

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
    int min_c = 7;
    int max_diff = 2;
    int min_span = 25;
    double vote_frac = 0.1;
    double vote_ratio = 3.0;
    bool help_requested = false;
    InputType input_type = InputType::AUTO;
};

struct Read {
    std::string seq_id;
    std::string seq;
    std::string qua;
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

            std::string _sid_buf;
            char _num[16];
            auto flush_sequence = [&]() {
                if (seq_id.empty() || seq.empty() || count >= depth) return;
                int seq_len = static_cast<int>(seq.size());
                for (int start = 0; start < seq_len && count < depth; start += step_size) {
                    int end = start + window_size;
                    if (end > seq_len) end = seq_len;
                    if (end - start < window_size) break;
                    std::string_view sv(seq.data() + start, end - start);
                    bool has_upper = false;
                    for (char c : sv) {
                        if (c >= 'A' && c <= 'Z') { has_upper = true; break; }
                    }
                    if (!has_upper) continue;
                    _sid_buf.clear();
                    _sid_buf += seq_id;
                    _sid_buf += ':';
                    auto rc = std::to_chars(_num, _num + 15, start + 1);
                    _sid_buf.append(_num, rc.ptr - _num);
                    _sid_buf += '-';
                    rc = std::to_chars(_num, _num + 15, end);
                    _sid_buf.append(_num, rc.ptr - _num);
                    Read read;
                    read.seq_id = _sid_buf;
                    read.seq = std::string(sv);
                    for (char& c : read.seq) {
                        if (c >= 'a' && c <= 'z') c = c - 'a' + 'A';
                    }
                    read.qua.clear();
                    chunk.push_back(std::move(read));
                    count++;
                    if (static_cast<int>(chunk.size()) >= chunk_size) {
                        flush_chunk(chunk);
                    }
                }
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

// ─── Worker thread ───

void worker_thread(int thread_id, ChunkQueue& queue, const MerQueryManager& mqm,
                   int mer_size, int min_c, int max_diff, int min_span,
                   double vote_frac, double vote_ratio,
                   const std::string& output_dir,
                   std::atomic<int64_t>& processed) {
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

    auto flush_buffer = [&]() {
        if (out_pos > 0) {
            fwrite(out_buf.data(), 1, out_pos, fp);
            out_pos = 0;
        }
    };

    std::vector<Read> chunk;
    while (queue.pop(chunk)) {
        for (const auto& read : chunk) {
            auto [confi_id, strand] = validate_read(read.seq, mqm, mer_size, min_c, max_diff, min_span, -1, vote_frac, vote_ratio);
            if (confi_id > -1) {
                int n = snprintf(out_buf.data() + out_pos, OUT_BUF_SIZE - out_pos,
                                 "%s\t%d\t%d\t%s\n",
                                 read.seq_id.c_str(), strand, confi_id, read.seq.c_str());
                if (n > 0) {
                    out_pos += static_cast<size_t>(n);
                }
                if (out_pos > OUT_BUF_SIZE * 0.8) {
                    flush_buffer();
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
        } else if (arg == "--min-c") {
            if (++i >= argc) throw std::runtime_error("Missing value for --min-c");
            args.min_c = std::stoi(argv[i]);
        } else if (arg == "--max-diff") {
            if (++i >= argc) throw std::runtime_error("Missing value for --max-diff");
            args.max_diff = std::stoi(argv[i]);
        } else if (arg == "--min-span") {
            if (++i >= argc) throw std::runtime_error("Missing value for --min-span");
            args.min_span = std::stoi(argv[i]);
        } else if (arg == "--vote-frac") {
            if (++i >= argc) throw std::runtime_error("Missing value for --vote-frac");
            args.vote_frac = std::stod(argv[i]);
        } else if (arg == "--vote-ratio") {
            if (++i >= argc) throw std::runtime_error("Missing value for --vote-ratio");
            args.vote_ratio = std::stod(argv[i]);
        } else if (arg[0] != '-') {
            args.reads_files.push_back(arg);
        } else {
            throw std::runtime_error("Unknown option: " + arg);
        }
        ++i;
    }
    return args;
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
                      << "  --min-c <n>           min passing adjacent pairs (default: 7)\n"
                      << "  --max-diff <n>        max gap diff tolerance (default: 2)\n"
                      << "  --min-span <n>        min read span covered (default: 25)\n"
                      << "  --vote-frac <f>       min fraction of k-mers for top CNE vote (default: 0.1)\n"
                      << "  --vote-ratio <r>      min ratio of 1st/2nd CNE votes (default: 3.0)\n";
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

        // Load mers table
        std::cerr << "Loading mers table from " << args.mers_file << " ...\n";
        MerQueryManager mqm;
        mqm.load_from_file(args.mers_file);
        int mer_size = mqm.get_mer_size();
        std::cerr << "  Loaded " << mqm.size() << " mers, k=" << mer_size << "\n";

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
        ChunkQueue queue;
        std::atomic<int64_t> processed{0};
        std::atomic<int64_t> total_reads{0};
        std::vector<std::thread> workers;

        for (int i = 0; i < args.threads; ++i) {
            workers.emplace_back(worker_thread, i, std::ref(queue),
                                 std::cref(mqm), mer_size,
                                 args.min_c, args.max_diff, args.min_span,
                                 args.vote_frac, args.vote_ratio,
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

        std::cerr << "Done. Output written to " << args.output_dir << "\n";
        return 0;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
