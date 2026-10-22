# CNEX Benchmark

多物种 genome benchmark，以及 reads 深度 / MSA 物种数消融实验。

## Active 目录

```text
benchmark/
  README.md
  Makefile
  benchmark.sh
  01_extract.py          # 采样 MSA，按 target 生成 queries/truth/blocks
  02_m8extract.py        # 每 block 取最佳命中
  03_analyze.py          # 多物种评估
  04_plot.py             # genome benchmark 图
  clean_queries.py       # query FASTA 清洗
  eval_locus.py          # reads locus 评估（Python）
  eval_locus.cpp         # eval_locus 源码（需编译）
  extract_m8_regions.py  # MMseqs m8 → regions
  reads_speed.py         # reads 计时采集
  run_reads_speed.sh     # CNEX reads validate+assemble 计时
  depth_experiment.py    # reads 深度主图 + S1/S2
  msa_ablation.py        # MSA 物种数消融
  msa_figure.py          # 附图 S3
  data/ dbs/ results/ stats/ figures/   # 生成物（gitignored）
```

## Genome benchmark

5 个 target：

```text
Homo_sapiens, Mus_musculus, Gallus_gallus,
Scyliorhinus_canicula, Callorhinchus_milii
```

```bash
make all
# 或分步
bash benchmark.sh extract Gallus_gallus
bash benchmark.sh dbs Gallus_gallus
bash benchmark.sh search Gallus_gallus
bash benchmark.sh m8extract Gallus_gallus
bash benchmark.sh analyze Gallus_gallus

python3 04_plot.py            # → figures/benchmark_multi.pdf
```

评估：命中与 truth 的 contig 名匹配，重叠 ≥30% 最短区间；MSA 用 seed=42
采样 10,000 blocks，CNEX mer 表 leave-one-out 构建（`-k 13 -c 4`）。

## Reads benchmark（head2head）

工作目录默认 `~/temp/head2head`，reads 为 0.25x–8x WGS 模拟数据。

```bash
CACHE=cold REPS=1 bash run_reads_speed.sh   # CNEX validate+assemble 计时
python3 reads_speed.py                      # → stats/reads_speed.tsv
python3 depth_experiment.py                 # → 主图 + S1/S2 + stats/depth_experiment.tsv
```

评估器 `eval_locus` 由 `eval_locus.cpp` 编译，依赖仓库的 `src/parasail`
submodule；二进制不进入 git。评估规则见 `stats/reads_locus_report.md`。

## MSA 物种数消融

人类 target，5x reads；4/8/12/16/19 个物种按 far → near 累积添加。

```bash
python3 msa_ablation.py       # 重建表 + validate + assemble（较慢）
python3 msa_figure.py         # → figures/msa_ablation_supp.pdf
```

4-species 点是 distant-only baseline；S3b 为 species composition 点阵。

## 图

```text
figures/benchmark_multi.pdf          # genome benchmark
figures/depth_experiment.pdf         # 主图：a recall / b recall vs cost
                                     #       c precision / d 5x memory
figures/depth_experiment_supp1.pdf   # S1：S-a–S-d
figures/depth_experiment_supp2.pdf   # S2：S-e–S-h
figures/msa_ablation_supp.pdf        # S3：MSA 物种数消融
```

## 软件与线程

- Megahit 1.2.9
- MMseqs2 18.8cc5c
- wgsim 1.0
- BLASTN 2.16.0
- SPAdes 4.3.0（已安装，未纳入结果）
- CNEX：仓库 commit + `bin/{validate,assemble,mertable}`

CNEX reads 计时用 `-t 4`，Megahit 用 `-t 64`。**core-hours 是主要计算量指标；
wall time 未按线程数归一化**，图和文中不要把它当作同线程对比。

## 数据与产物

`data/`、`dbs/`、`results/`、`stats/`、`figures/` 均为生成物并已 gitignore。
`src/parasail` 是 submodule；编译产物不提交。旧版 querylib / megahit_test /
backup 数据已移到 repo 外的
`/43T/home/zhengjiangmin/Software/CNEX_benchmark_archive_20261022/`。
