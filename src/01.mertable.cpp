#include <hip/mertable_core.hpp>
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <set>
#include <unordered_map>
#include <vector>

struct Args {
    std::string msa_file;
    int k = 13;
    int min_c = 4;
    float min_entropy = 1.2f;
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

// Split MSA by ### blocks
std::vector<std::string> read_bunches(const std::string& path) {
    std::ifstream f(path);
    if (!f.is_open()) throw std::runtime_error("Cannot open: " + path);
    std::string buf((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
    std::vector<std::string> bunches;
    size_t pos = 0;
    while (true) {
        auto sep = buf.find("###\n", pos);
        if (sep == std::string::npos) {
            auto remainder = buf.substr(pos);
            if (!remainder.empty()) {
                // trim
                while (!remainder.empty() && (remainder.back() == '\n' || remainder.back() == ' '))
                    remainder.pop_back();
                if (!remainder.empty())
                    bunches.push_back(remainder);
            }
            break;
        }
        auto block = buf.substr(pos, sep - pos);
        if (!block.empty()) {
            while (!block.empty() && block.back() == '\n') block.pop_back();
            if (!block.empty())
                bunches.push_back(block);
        }
        pos = sep + 4;
    }
    return bunches;
}

// Parse a bunch into sequence dict, skip "ref" sequences and gap-only sequences
std::unordered_map<std::string, std::string> bunch2fas(const std::string& bunch) {
    std::unordered_map<std::string, std::string> fas;
    std::istringstream ss(bunch);
    std::string line, seq_id, seq;
    while (std::getline(ss, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!seq_id.empty() && !seq.empty()) {
                bool has_base = false;
                for (char c : seq) if (c != '-') { has_base = true; break; }
                if (has_base) fas[seq_id] = seq;
            }
            seq_id = line.substr(1);
            // trim whitespace from seq_id
            while (!seq_id.empty() && (seq_id.back() == ' ' || seq_id.back() == '\t' || seq_id.back() == '\r'))
                seq_id.pop_back();
            seq.clear();
        } else {
            for (char c : line) {
                if (c != ' ' && c != '\r') seq.push_back(static_cast<char>(std::toupper(c)));
            }
        }
    }
    if (!seq_id.empty() && !seq.empty()) {
        bool has_base = false;
        for (char c : seq) if (c != '-') { has_base = true; break; }
        if (has_base) fas[seq_id] = seq;
    }
    return fas;
}

// Gap-aware k-mer generation (same as Python generate_mers)
struct MerPos { int loci; std::string mer; };
std::vector<MerPos> generate_mers(const std::string& seq, int mer_size) {
    std::vector<MerPos> result;
    std::string pure_seq;
    for (char c : seq) if (c != '-') pure_seq.push_back(c);
    if ((int)pure_seq.size() < mer_size * 2) return result;

    // Build gap bias map: pure_seq index -> seq index
    std::unordered_map<int, int> bias;
    int bias_val = 0;
    for (int i = 0; i < (int)seq.size(); ++i) {
        if (seq[i] == '-') {
            bias_val++;
        } else {
            bias[i - bias_val] = i;
        }
    }

    for (int i = 0; i + mer_size <= (int)pure_seq.size(); ++i) {
        auto local_seq = pure_seq.substr(i, mer_size);
        bool valid = true;
        for (char c : local_seq) if (c == ' ') { valid = false; break; }
        if (valid) {
            auto it = bias.find(i);
            if (it != bias.end())
                result.push_back({it->second, local_seq});
        }
    }
    return result;
}

int main(int argc, char* argv[]) {
    try {
        Args args = parse_args(argc, argv);

        if (args.help_requested || args.msa_file.empty()) {
            std::cerr << "Usage: " << argv[0]
                      << " <msa_file> [options]\n"
                      << "\nOptions:\n"
                      << "  -k, --mer-size <n>    K-mer size (default: 13)\n"
                      << "  -c, --min-c <n>       Minimum species count (default: 4)\n"
                      << "  --min-entropy <f>     Minimum k-mer entropy (default: 1.2)\n"
                      << "  -o, --output <file>    Output TSV (default: mers_table.tsv)\n";
            return 1;
        }

        if (args.k > 16) throw std::runtime_error("k must be <= 16");
        if (args.k < 7) throw std::runtime_error("k must be >= 7");

        std::cerr << "Reading MSA from " << args.msa_file << " ...\n";
        auto bunches = read_bunches(args.msa_file);
        std::cerr << "  Total bunches: " << bunches.size() << "\n";

        // Open ref.fa output
        std::ofstream ref_out;
        {
            auto pos = args.output.rfind('/');
            std::string ref_path = (pos != std::string::npos)
                ? args.output.substr(0, pos + 1) + "ref.fa" : "ref.fa";
            ref_out.open(ref_path);
            if (ref_out.is_open())
                std::cerr << "  Writing ref sequences to " << ref_path << "\n";
        }

        TableMer mertable(args.k);
        mertable.set_min_entropy(args.min_entropy);

        for (size_t bunch_id = 0; bunch_id < bunches.size(); ++bunch_id) {
            if (bunch_id % 1000 == 0)
                std::cerr << "\r  Processing bunch " << (bunch_id + 1) << "/" << bunches.size();
            auto fas = bunch2fas(bunches[bunch_id]);
            // Count real species (excluding ref)
            int real_species = 0;
            for (const auto& [header, _] : fas) {
                if (header.size() < 3 || header.substr(0, 3) != "ref") real_species++;
            }
            if (real_species < args.min_c) continue;

            // Write ref sequence: median-length non-ref species, de-gapped
            if (ref_out.is_open()) {
                std::vector<std::pair<size_t, std::string>> cand;
                for (const auto& [hdr, seq] : fas) {
                    if (hdr.size() >= 3 && hdr.substr(0, 3) == "ref") continue;
                    std::string raw;
                    for (char c : seq) if (c != '-') raw.push_back(c);
                    cand.emplace_back(raw.size(), std::move(raw));
                }
                if (!cand.empty()) {
                    std::sort(cand.begin(), cand.end());
                    ref_out << ">" << bunch_id << "\n" << cand[cand.size() * 3 / 4].second << "\n";
                }
            }

            // Track (mer -> (bunch_id, loci)) pairs
            std::unordered_map<std::string, std::vector<std::pair<int, int>>> local_affi;
            std::unordered_map<std::string, int> local_mer_count;

            for (const auto& [header, seq] : fas) {
                if (header.size() >= 3 && header.substr(0, 3) == "ref") continue;
                for (const auto& mp : generate_mers(seq, args.k)) {
                    local_affi[mp.mer].push_back({static_cast<int>(bunch_id), mp.loci});
                    local_mer_count[mp.mer]++;
                }
            }

            for (const auto& [mer, locations] : local_affi) {
                // Check if this mer maps to unique (bunch_id, loci) pair
                std::set<std::pair<int, int>> unique_locs(locations.begin(), locations.end());
                if (unique_locs.size() > 1) continue;
                auto [bid, loci] = *unique_locs.begin();
                int count = local_mer_count[mer];
                mertable.add(mer, bid, loci, count);
            }
        }

        std::cerr << "\n  Writing " << mertable.size() << " mers to " << args.output << " ...\n";
        mertable.dump(args.output);
        std::cerr << "Done.\n";
        return 0;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
