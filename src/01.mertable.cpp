#include <hip/mertable_core.hpp>
#include <algorithm>
#include <array>
#include <atomic>
#include <charconv>
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <string_view>
#include <thread>
#include <unordered_map>
#include <vector>

struct Args {
    std::string msa_file;
    int k = 13;
    int min_c = 4;
    float min_entropy = 1.2f;
    int max_assign = 0;   // 0 = keep all assignments
    int threads = 0;      // 0 = auto
    std::string output = "mers_table.tsv";
    bool help_requested = false;
};

Args parse_args(int argc, char* argv[]) {
    Args args;
    int i = 1;
    while (i < argc) {
        std::string arg(argv[i]);
        if (arg == "-h" || arg == "--help") {
            args.help_requested = true;
            return args;
        } else if (arg == "-k" || arg == "--mer-size") {
            if (++i >= argc) throw std::runtime_error("Missing value for -k");
            args.k = std::stoi(argv[i]);
        } else if (arg == "-c" || arg == "--min-c") {
            if (++i >= argc) throw std::runtime_error("Missing value for -c");
            args.min_c = std::stoi(argv[i]);
        } else if (arg == "--min-entropy") {
            if (++i >= argc) throw std::runtime_error("Missing value for --min-entropy");
            args.min_entropy = std::stof(argv[i]);
        } else if (arg == "-a" || arg == "--max-assign") {
            if (++i >= argc) throw std::runtime_error("Missing value for -a");
            args.max_assign = std::stoi(argv[i]);
        } else if (arg == "-t" || arg == "--threads") {
            if (++i >= argc) throw std::runtime_error("Missing value for -t");
            args.threads = std::stoi(argv[i]);
        } else if (arg == "-o" || arg == "--output") {
            if (++i >= argc) throw std::runtime_error("Missing value for -o");
            args.output = argv[i];
        } else if (arg[0] != '-') {
            if (args.msa_file.empty()) args.msa_file = arg;
            else throw std::runtime_error("Unexpected argument: " + arg);
        } else {
            throw std::runtime_error("Unknown option: " + arg);
        }
        ++i;
    }
    return args;
}

static std::string read_entire_file(const std::string& path) {
    FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) throw std::runtime_error("Cannot open: " + path);
    if (std::fseek(f, 0, SEEK_END) != 0) {
        std::fclose(f);
        throw std::runtime_error("Cannot seek: " + path);
    }
    long sz = std::ftell(f);
    if (sz < 0) {
        std::fclose(f);
        throw std::runtime_error("Cannot stat size: " + path);
    }
    std::rewind(f);
    std::string buf(static_cast<size_t>(sz), '\0');
    size_t got = 0;
    while (got < buf.size()) {
        size_t n = std::fread(&buf[got], 1, buf.size() - got, f);
        if (n == 0) break;
        got += n;
    }
    std::fclose(f);
    buf.resize(got);
    return buf;
}

struct BunchSlice {
    size_t off;
    size_t len;
};

static std::vector<BunchSlice> slice_bunches(const std::string& buf) {
    std::vector<BunchSlice> bunches;
    bunches.reserve(buf.size() / 1024 + 16);
    size_t pos = 0;
    while (true) {
        size_t sep = buf.find("###\n", pos);
        if (sep == std::string::npos) {
            size_t len = buf.size() - pos;
            while (len > 0 && (buf[pos + len - 1] == '\n' || buf[pos + len - 1] == ' ')) {
                --len;
            }
            if (len > 0) bunches.push_back({pos, len});
            break;
        }
        size_t len = sep - pos;
        while (len > 0 && buf[pos + len - 1] == '\n') {
            --len;
        }
        if (len > 0) bunches.push_back({pos, len});
        pos = sep + 4;
    }
    return bunches;
}

struct MerEntry {
    uint32_t code;
    int bunch_id;
    int loci;
    int count;
};

static constexpr size_t kNumShards = 64;

int main(int argc, char* argv[]) {
    try {
        Args args = parse_args(argc, argv);

        if (args.help_requested || args.msa_file.empty()) {
            std::cerr << "Usage: " << argv[0]
                      << " <msa_file> [options]\n"
                      << "\nOptions:\n"
                      << "  -k, --mer-size <n>    K-mer size (default: 13)\n"
                      << "  -c, --min-c <n>       Minimum species count (default: 4)\n"
                      << "  -a, --max-assign <n>  Optional build-time cap: drop k-mers with more than n\n"
                      << "                        (element, locus) assignments (default: 0 = keep all;\n"
                      << "                        filtering is normally done at validate load time)\n"
                      << "  --min-entropy <f>     Minimum k-mer entropy (default: 1.2)\n"
                      << "  -t, --threads <n>     Worker threads (default: auto up to 16)\n"
                      << "  -o, --output <file>   Output TSV (default: mers_table.tsv)\n";
            return 1;
        }

        if (args.k > 16) throw std::runtime_error("k must be <= 16");
        if (args.k < 7) throw std::runtime_error("k must be >= 7");

        int nthreads = args.threads;
        if (nthreads <= 0) {
            unsigned hw = std::thread::hardware_concurrency();
            nthreads = static_cast<int>(std::min<unsigned>(16u, std::max<unsigned>(1u, hw)));
        }

        std::cerr << "Reading MSA from " << args.msa_file << " ...\n";
        std::string buf = read_entire_file(args.msa_file);
        auto bunches = slice_bunches(buf);
        std::cerr << "  Total bunches: " << bunches.size() << "\n";

        const int k = args.k;
        const int min_c = args.min_c;
        const float min_entropy = args.min_entropy;
        const uint32_t mer_mask = static_cast<uint32_t>((1ULL << (2 * k)) - 1ULL);
        const int shard_shift = 2 * k - 6;

        static const std::array<int8_t, 256> base_lut = [] {
            std::array<int8_t, 256> t{};
            t.fill(-1);
            t['A'] = 0; t['C'] = 1; t['G'] = 2; t['T'] = 3;
            t['a'] = 0; t['c'] = 1; t['g'] = 2; t['t'] = 3;
            return t;
        }();

        std::vector<std::array<std::vector<MerEntry>, kNumShards>> tls_shards(nthreads);
        std::atomic<size_t> next_bunch{0};
        const size_t n_bunches = bunches.size();
        const size_t block_step = 256;

        auto worker = [&](int tid) {
            auto& my_shards = tls_shards[tid];
            struct SeqSpan {
                std::string_view id;
                size_t off;
                size_t len;
            };
            std::string seq_pool;
            std::vector<SeqSpan> seqs;
            std::vector<SeqSpan> dedup_seqs;
            std::vector<char> pure;
            std::vector<int> pure_loci;
            std::vector<std::pair<uint32_t, int>> bunch_mers;
            seq_pool.reserve(65536);
            seqs.reserve(64);
            dedup_seqs.reserve(64);
            pure.reserve(4096);
            pure_loci.reserve(4096);
            bunch_mers.reserve(8192);

            while (true) {
                size_t b0 = next_bunch.fetch_add(block_step, std::memory_order_relaxed);
                if (b0 >= n_bunches) break;
                size_t b1 = std::min(n_bunches, b0 + block_step);

                for (size_t bunch_id = b0; bunch_id < b1; ++bunch_id) {
                    const char* p = buf.data() + bunches[bunch_id].off;
                    const char* end = p + bunches[bunch_id].len;

                    seq_pool.clear();
                    seqs.clear();

                    std::string_view cur_id;
                    size_t cur_off = 0;
                    bool cur_has_base = false;

                    auto flush_seq = [&]() {
                        if (!cur_id.empty() && seq_pool.size() > cur_off && cur_has_base) {
                            seqs.push_back({cur_id, cur_off, seq_pool.size() - cur_off});
                        }
                    };

                    while (p < end) {
                        const void* nl_ptr = std::memchr(p, '\n', static_cast<size_t>(end - p));
                        const char* line_end = nl_ptr ? static_cast<const char*>(nl_ptr) : end;
                        size_t line_len = static_cast<size_t>(line_end - p);
                        if (line_len > 0) {
                            if (p[0] == '>') {
                                flush_seq();
                                const char* id_beg = p + 1;
                                const char* id_end = line_end;
                                while (id_end > id_beg &&
                                       (id_end[-1] == ' ' || id_end[-1] == '\t' || id_end[-1] == '\r')) {
                                    --id_end;
                                }
                                cur_id = std::string_view(id_beg, static_cast<size_t>(id_end - id_beg));
                                cur_off = seq_pool.size();
                                cur_has_base = false;
                            } else {
                                for (const char* q = p; q < line_end; ++q) {
                                    unsigned char uc = static_cast<unsigned char>(*q);
                                    if (uc == ' ' || uc == '\r') continue;
                                    char c = static_cast<char>(std::toupper(uc));
                                    if (c != '-') cur_has_base = true;
                                    seq_pool.push_back(c);
                                }
                            }
                        }
                        p = nl_ptr ? (line_end + 1) : end;
                    }
                    flush_seq();

                    if (seqs.empty()) continue;

                    // Deduplicate duplicate headers in the same bunch by keeping the last
                    // non-gap-only sequence, matching unordered_map overwrite semantics.
                    dedup_seqs.clear();
                    if (seqs.size() <= 32) {
                        for (size_t si = 0; si < seqs.size(); ++si) {
                            bool overwritten = false;
                            for (size_t sj = si + 1; sj < seqs.size(); ++sj) {
                                if (seqs[si].id == seqs[sj].id) {
                                    overwritten = true;
                                    break;
                                }
                            }
                            if (!overwritten) dedup_seqs.push_back(seqs[si]);
                        }
                    } else {
                        std::unordered_map<std::string_view, size_t> last_idx;
                        last_idx.reserve(seqs.size());
                        for (size_t si = 0; si < seqs.size(); ++si) {
                            last_idx[seqs[si].id] = si;
                        }
                        for (size_t si = 0; si < seqs.size(); ++si) {
                            if (last_idx[seqs[si].id] == si) dedup_seqs.push_back(seqs[si]);
                        }
                    }

                    int real_species = 0;
                    for (const auto& sp : dedup_seqs) {
                        if (sp.id.size() < 3 || sp.id.compare(0, 3, "ref") != 0) {
                            ++real_species;
                        }
                    }
                    if (real_species < min_c) continue;

                    bunch_mers.clear();
                    for (const auto& sp : dedup_seqs) {
                        if (sp.id.size() >= 3 && sp.id.compare(0, 3, "ref") == 0) continue;
                        pure.clear();
                        pure_loci.clear();
                        const char* s = seq_pool.data() + sp.off;
                        for (size_t idx = 0; idx < sp.len; ++idx) {
                            char c = s[idx];
                            if (c != '-') {
                                pure.push_back(c);
                                pure_loci.push_back(static_cast<int>(idx));
                            }
                        }
                        if (static_cast<int>(pure.size()) < k * 2) continue;

                        uint32_t rolling = 0;
                        int valid_run = 0;
                        const int pure_n = static_cast<int>(pure.size());
                        for (int idx = 0; idx < pure_n; ++idx) {
                            int8_t bits = base_lut[static_cast<unsigned char>(pure[idx])];
                            if (bits < 0) {
                                if (min_entropy <= 0.0f) {
                                    // Legacy fallback for invalid base when min_entropy <= 0: mer2int returns 0
                                    if (idx + 1 >= k) {
                                        bunch_mers.emplace_back(0u, pure_loci[idx - k + 1]);
                                    }
                                }
                                valid_run = 0;
                                rolling = 0;
                                continue;
                            }
                            rolling = ((rolling << 2) | static_cast<uint32_t>(bits)) & mer_mask;
                            ++valid_run;
                            if (valid_run >= k) {
                                bunch_mers.emplace_back(rolling, pure_loci[idx - k + 1]);
                            }
                        }
                    }

                    if (bunch_mers.empty()) continue;
                    std::sort(bunch_mers.begin(), bunch_mers.end(),
                              [](const std::pair<uint32_t, int>& a, const std::pair<uint32_t, int>& b) {
                                  return a.first < b.first;
                              });

                    const size_t bm_n = bunch_mers.size();
                    size_t r0 = 0;
                    while (r0 < bm_n) {
                        uint32_t code = bunch_mers[r0].first;
                        int loci0 = bunch_mers[r0].second;
                        size_t r1 = r0 + 1;
                        bool unique_loc = true;
                        while (r1 < bm_n && bunch_mers[r1].first == code) {
                            if (bunch_mers[r1].second != loci0) unique_loc = false;
                            ++r1;
                        }
                        if (unique_loc && mer_entropy_code(code, k) >= min_entropy) {
                            int cnt = static_cast<int>(r1 - r0);
                            uint32_t shard = (code >> shard_shift) & (kNumShards - 1);
                            my_shards[shard].push_back(
                                MerEntry{code, static_cast<int>(bunch_id), loci0, cnt});
                        }
                        r0 = r1;
                    }
                }
            }
        };

        std::vector<std::thread> workers;
        workers.reserve(nthreads);
        for (int t = 0; t < nthreads; ++t) workers.emplace_back(worker, t);
        for (auto& w : workers) w.join();
        // Free raw MSA text buffer before shard merge to keep peak RSS low.
        {
            std::string().swap(buf);
            std::vector<BunchSlice>().swap(bunches);
        }

        // Parallel merge & TSV formatting per shard
        std::array<std::string, kNumShards> shard_text;
        std::atomic<size_t> next_shard{0};
        std::atomic<size_t> total_mers{0};

        auto merge_worker = [&]() {
            static const char INT_TO_BASE[4] = {'A', 'C', 'G', 'T'};
            std::vector<MerEntry> merged;
            char mer_buf[16];
            char num_buf[64];

            while (true) {
                size_t s = next_shard.fetch_add(1, std::memory_order_relaxed);
                if (s >= kNumShards) break;

                size_t total_entries = 0;
                for (int t = 0; t < nthreads; ++t) total_entries += tls_shards[t][s].size();
                if (total_entries == 0) continue;

                merged.clear();
                merged.reserve(total_entries);
                for (int t = 0; t < nthreads; ++t) {
                    auto& v = tls_shards[t][s];
                    merged.insert(merged.end(), v.begin(), v.end());
                    std::vector<MerEntry>().swap(v);
                }

                std::sort(merged.begin(), merged.end(), [](const MerEntry& a, const MerEntry& b) {
                    if (a.code != b.code) return a.code < b.code;
                    if (a.bunch_id != b.bunch_id) return a.bunch_id < b.bunch_id;
                    return a.loci < b.loci;
                });

                std::string& out = shard_text[s];
                out.reserve(merged.size() * 28);
                size_t local_mers = 0;

                size_t g0 = 0;
                const size_t mn = merged.size();
                while (g0 < mn) {
                    uint32_t code = merged[g0].code;
                    size_t g1 = g0 + 1;
                    while (g1 < mn && merged[g1].code == code) ++g1;

                    if (args.max_assign == 0 || static_cast<int>(g1 - g0) <= args.max_assign) {
                        ++local_mers;
                        uint32_t tmp = code;
                        for (int b = k - 1; b >= 0; --b) {
                            mer_buf[b] = INT_TO_BASE[tmp & 3u];
                            tmp >>= 2;
                        }
                        for (size_t idx = g0; idx < g1; ++idx) {
                            out.append(mer_buf, static_cast<size_t>(k));
                            out.push_back('\t');
                            auto r1 = std::to_chars(num_buf, num_buf + 20, merged[idx].bunch_id);
                            out.append(num_buf, static_cast<size_t>(r1.ptr - num_buf));
                            out.push_back('\t');
                            auto r2 = std::to_chars(num_buf, num_buf + 20, merged[idx].loci);
                            out.append(num_buf, static_cast<size_t>(r2.ptr - num_buf));
                            out.push_back('\t');
                            auto r3 = std::to_chars(num_buf, num_buf + 20, merged[idx].count);
                            out.append(num_buf, static_cast<size_t>(r3.ptr - num_buf));
                            out.push_back('\n');
                        }
                    }
                    g0 = g1;
                }
                total_mers.fetch_add(local_mers, std::memory_order_relaxed);
            }
        };

        workers.clear();
        for (int t = 0; t < nthreads; ++t) workers.emplace_back(merge_worker);
        for (auto& w : workers) w.join();

        std::cerr << "\n  Writing " << total_mers.load() << " mers to " << args.output << " ...\n";
        FILE* fout = std::fopen(args.output.c_str(), "wb");
        if (!fout) throw std::runtime_error("Cannot open output file: " + args.output);
        static char write_buf[1 << 20];
        std::setvbuf(fout, write_buf, _IOFBF, sizeof(write_buf));
        for (size_t s = 0; s < kNumShards; ++s) {
            if (!shard_text[s].empty()) {
                if (std::fwrite(shard_text[s].data(), 1, shard_text[s].size(), fout) != shard_text[s].size()) {
                    std::fclose(fout);
                    throw std::runtime_error("Write failed for " + args.output);
                }
                std::string().swap(shard_text[s]);
            }
        }
        std::fclose(fout);
        std::cerr << "Done.\n";
        return 0;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
