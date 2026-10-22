// Locus-referenced evaluation (C++ port of benchmark/eval_locus.py).
//
// A claim (CNEX contig `>ele` or megahit region `>block_XXXXX|...`) is TP iff
// a parasail SW alignment (forward or reverse-complement) against any
// truth-locus reference sequence of the SAME block has identity >= THR and
// aligned length >= COV * min(len(claim), len(ref)).
//
// Usage:
//   eval_locus <refs.fa> <claims.fa> [thr=0.95] [cov=0.30]
//
// refs.fa: headers are block ids (integers), as written by eval_locus.py
//          build-refs. claims.fa: FASTA.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>

#include "parasail.h"

static const char* COMP = "TGCAN";  // ACGTN -> TGCAN

static std::string revcomp(const std::string& s) {
    std::string r(s.size(), 'N');
    for (size_t i = 0; i < s.size(); ++i) {
        char c = s[s.size() - 1 - i];
        switch (c) {
            case 'A': r[i] = 'T'; break; case 'C': r[i] = 'G'; break;
            case 'G': r[i] = 'C'; break; case 'T': r[i] = 'A'; break;
            case 'N': r[i] = 'N'; break; default: r[i] = 'N'; break;
        }
    }
    return r;
}

static void load_fasta(const char* path, std::vector<std::pair<int, std::string>>& out) {
    FILE* f = fopen(path, "r");
    if (!f) { fprintf(stderr, "cannot open %s\n", path); exit(2); }
    char* line = nullptr; size_t cap = 0; ssize_t len;
    int cur = -1; std::string seq;
    auto flush = [&]() { if (cur >= 0 && !seq.empty()) out.emplace_back(cur, seq); seq.clear(); };
    while ((len = getline(&line, &cap, f)) != -1) {
        if (len > 0 && line[len - 1] == '\n') line[--len] = 0;
        if (len == 0) continue;
        if (line[0] == '>') {
            flush();
            const char* h = line + 1;
            if (strncmp(h, "block_", 6) == 0) cur = atoi(h + 6);
            else cur = atoi(h);  // leading integer up to '.'/space
        } else {
            for (ssize_t i = 0; i < len; ++i) {
                char c = line[i];
                if (c != ' ' && c != '\r') seq.push_back(c);
            }
        }
    }
    flush();
    free(line);
    fclose(f);
}

int main(int argc, char** argv) {
    if (argc < 3) {
        fprintf(stderr, "usage: %s <refs.fa> <claims.fa> [thr=0.95] [cov=0.30]\n", argv[0]);
        return 1;
    }
    double thr = (argc > 3) ? atof(argv[3]) : 0.95;
    double cov = (argc > 4) ? atof(argv[4]) : 0.30;

    std::vector<std::pair<int, std::string>> refv, claimv;
    load_fasta(argv[1], refv);
    load_fasta(argv[2], claimv);

    std::unordered_map<int, std::vector<const std::string*>> refs;
    for (auto& p : refv) refs[p.first].push_back(&p.second);
    int nt = (int)refs.size();

    parasail_matrix_t* mat = parasail_matrix_create("ACGTN", 2, -1);

    long long tp = 0, fp = 0, in_truth = 0;
    std::unordered_map<int, char> blocks;
    struct AuditRow { size_t claim_idx; int ele, qlen, ok, ref_idx, strand; double ident; int alen; };
    std::vector<AuditRow> audit;
    size_t claim_idx = 0;
    for (auto& c : claimv) {
        const size_t this_claim = claim_idx++;
        int ele = c.first;
        auto it = refs.find(ele);
        if (it == refs.end()) continue;
        ++in_truth;
        const std::string& qf = c.second;
        std::string qr = revcomp(qf);
        bool ok = false;
        double best_ident = 0.0;
        int best_alen = 0, best_ref = -1, best_strand = 0;
        for (size_t ri = 0; ri < it->second.size(); ++ri) {
            const std::string& rs = *it->second[ri];
            int rslen = (int)rs.size();
            const std::string* qs[2] = { &qf, &qr };
            for (int si = 0; si < 2; ++si) {
                const std::string* q = qs[si];
                int qlen = (int)q->size();
                if (qlen == 0 || rslen == 0) continue;
                parasail_result_t* res = parasail_sw_stats_striped_16(
                    q->c_str(), qlen, rs.c_str(), rslen, 1, 1, mat);
                if (res->stats->length > 0) {
                    double ident = (double)res->stats->matches / (double)res->stats->length;
                    int alen = (int)res->stats->length;
                    if (ident > best_ident || (ident == best_ident && alen > best_alen)) {
                        best_ident = ident; best_alen = alen;
                        best_ref = (int)ri; best_strand = si ? -1 : 1;
                    }
                    int mn = std::min(qlen, rslen);
                    if (ident >= thr && res->stats->length >= cov * mn) ok = true;
                }
                parasail_result_free(res);
                if (ok) break;
            }
            if (ok) break;
        }
        tp += ok;
        fp += !ok;
        blocks[ele] = ok ? 1 : 0;
        audit.push_back({this_claim, ele, (int)qf.size(), (int)ok, best_ref, best_strand, best_ident, best_alen});
    }
    long long blocks_tp = 0;
    for (auto& kv : blocks) blocks_tp += kv.second;

    double recall = nt ? (double)blocks_tp / nt * 100.0 : 0.0;
    double prec = (tp + fp) ? (double)tp / (tp + fp) * 100.0 : 0.0;
    printf("claims=%zu in_truth=%lld tp=%lld fp=%lld blocks_tp=%lld nt=%d "
           "recall=%.1f%% prec=%.1f%% (thr=%.2f cov=%.2f)\n",
           claimv.size(), in_truth, tp, fp, blocks_tp, nt, recall, prec, thr, cov);

    if (argc > 5) {
        FILE* df = fopen(argv[5], "w");
        if (df) {
            for (auto& c : claimv) {
                if (c.first < 0) continue;
                fprintf(df, "%d\t%d\n", c.first, (int)blocks[c.first]);
            }
            fclose(df);
        }
    }

    if (argc > 6) {
        FILE* af = fopen(argv[6], "w");
        if (af) {
            fprintf(af, "claim_idx\tele\tqlen\tok\tbest_ident\tbest_alen\tbest_ref_idx\tbest_strand\n");
            for (auto& r : audit)
                fprintf(af, "%zu\t%d\t%d\t%d\t%.4f\t%d\t%d\t%d\n",
                        r.claim_idx, r.ele, r.qlen, r.ok, r.ident, r.alen, r.ref_idx, r.strand);
            fclose(af);
        }
    }

    parasail_matrix_free(mat);
    return 0;
}
