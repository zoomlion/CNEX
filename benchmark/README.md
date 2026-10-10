# CNEX Benchmark

多物种 benchmark:CNEX vs MMseqs2 vs BLASTN,5 个目标基因组
(Homo_sapiens, Mus_musculus, Gallus_gallus, Xenopus_tropicalis, Callorhinchus_milii)。

## 原理

从 `most-cons-cne.filtered.fa`(26k 块多物种 MSA)按 seed=42 采样 10,000 块。
对每个目标物种 T:
- queries = 所有非 T、非 ref 物种的序列(149k 量级)
- ground truth = T 在各块的基因组坐标(过滤:非-gap ≥ 25bp)
- CNEX mer 表由"剔除 T"的比对块构建(leave-one-out,公平的同源搜索)
- 评估:命中与 truth 的 **contig 名精确匹配** + 重叠 ≥30% 最短区间

## Usage

```bash
make all                # 完整流程(5 物种)
bash benchmark.sh all   # 同上

# 分步(可指定物种子集)
bash benchmark.sh extract Gallus_gallus        # 采样+按物种出 queries/truth/blocks
bash benchmark.sh dbs Gallus_gallus            # 建 3 个库
bash benchmark.sh search Gallus_gallus         # CNEX + MMseqs2 + BLASTN
bash benchmark.sh m8extract Gallus_gallus      # 每 block 最佳命中
bash benchmark.sh sensitivity                  # MMseqs2 s=1.0,4.0,7.5,5.7(仅人类)
bash benchmark.sh querylib Homo_sapiens        # CNEX query-库子采样 sweep (10/20/40/60/80%+100%)
bash benchmark.sh analyze Gallus_gallus        # 汇总 results → stats/summary.tsv
python3 04_plot.py                             # figures/benchmark_multi.pdf

make clean              # 清理衍生数据(backup_human_* 保留)
```

querylib 分步: `querylib-extract` / `querylib-search` / `querylib-analyze` / `querylib-plot`

## 目录结构

```
data/genomes/{Species}/{Species}.fa(.fai)   # 目标基因组(RefSeq 命名,与 MSA 版本一致)
data/queries/{Species}.fa                   # 物种专属 query
data/truth/{Species}.tsv                    # 物种专属 ground truth
data/blocks/{Species}.fa                    # 剔除目标物种的比对块(CNEX mer 表输入)
data/querylib/{Species}/bf{frac}/           # query-库子采样产物(queries/blocks/kept_blocks.tsv)
dbs/{Species}/{mmseqs,blast,cnex}/          # 每物种数据库
results/{Species}/{mmseqs,blastn,cnex}/     # 每物种搜索结果
results/querylib/{Species}/bf{frac}/cnex/   # 子采样级别搜索结果(CNEX)
stats/{Species}/time_*.txt                  # /usr/bin/time -v 日志
stats/summary.tsv                           # 汇总表(每物种 × 3 工具)
stats/querylib_summary.tsv                  # query-库 sweep 汇总
figures/benchmark_multi.pdf                 # 统计图
figures/querylib_recall.pdf                 # query-库 sweep 图
backup_human_20260801/                      # 旧单物种(人类)产物备份
```

## 文件说明

- `01_extract.py` — 采样、按 `--target` 提取 queries/truth/blocks(支持 `--block-frac` 嵌套子采样)
- `02_m8extract.py` — 每物种、每 block 取 evalue<1e-5 最佳命中
- `03_analyze.py` — 多物种评估、出表
- `04_plot.py` — 多物种统计图
- `05_querylib_analyze.py` — query-库子采样评估(相对/绝对 recall、恢复率、bootstrap CI)
- `06_querylib_plot.py` — query-库子采样图
- `benchmark.sh` — 统一调用入口
