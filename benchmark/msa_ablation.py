#!/usr/bin/env python3
"""MSA species-count ablation at 5x (human target).

For each species count N in 4/8/12/16/19, subset the leave-one-out human
block MSA to a cumulative far-to-near species set, rebuild the
13-mer table with the production min_c=4, then run validate + assemble at 5x
and evaluate with eval_locus.

Usage:
    python3 benchmark/msa_ablation.py [--force]
Outputs:
    /tmp/zjm_cnex8/msa_ablation/           intermediate files
    benchmark/stats/msa_ablation.tsv       metrics
"""
import argparse
import os
import re
import shutil
import subprocess
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
BASE = Path(os.environ.get('BASE', Path.home() / 'temp' / 'head2head'))
WORK = Path('/tmp/zjm_cnex8/msa_ablation')
MSA = HERE / 'data' / 'blocks' / 'Homo_sapiens.fa'
MERTABLE = ROOT / 'mertable'
VALIDATE = ROOT / 'validate'
ASSEMBLE = ROOT / 'assemble'
PIGZ = BASE / 'bin' / 'pigz'
REFS = BASE / 'eval' / 'locus' / 'refs.fa'
EVAL = HERE / 'eval_locus'
OUT = HERE / 'stats' / 'msa_ablation.tsv'

# Cumulative, far-to-near order (most distant -> closest relative of human).
# Human is the leave-one-out target and is absent from the block MSA.
SPECIES = [
    'Callorhinchus_milii', 'Rhincodon_typus', 'Chiloscyllium_plagiosum',
    'Carcharodon_carcharias', 'Pristis_pectinata', 'Mobula_birostris',
    'Hypanus_sabinus', 'Leucoraja_erinacea', 'Amblyraja_radiata',
    'Scyliorhinus_canicula', 'Polypterus_senegalus', 'Danio_rerio',
    'Latimeria_chalumnae', 'Bombina_bombina', 'Xenopus_tropicalis',
    'Gallus_gallus', 'Gracilinanus_agilis', 'Dromiciops_gliroides',
    'Mus_musculus',
]
LEVELS = [4, 8, 12, 16, 19]


def species_of(header):
    return header[1:].split('.')[0]


def subset_blocks(level, out_path):
    """Write a block MSA containing only the first `level` species."""
    keep = set(SPECIES[:level])
    blocks = []
    cur = []
    for line in open(MSA):
        if line.rstrip('\n') == '###':
            if cur:
                blocks.append(cur)
            cur = []
        else:
            cur.append(line)
    if cur:
        blocks.append(cur)

    with open(out_path, 'w') as f:
        written = 0
        for block in blocks:
            filtered = []
            i = 0
            while i < len(block):
                line = block[i]
                if line.startswith('>'):
                    hdr = line
                    seq = block[i + 1] if i + 1 < len(block) else ''
                    if species_of(hdr) in keep:
                        filtered.append(hdr)
                        filtered.append(seq)
                    i += 2
                else:
                    i += 1
            # Preserve the block index: the table builder numbers blocks
            # sequentially, so an omitted block would shift every element id.
            if not filtered:
                filtered = ['>ref_dummy\n', 'N\n']
            if written:
                f.write('###\n')
            f.writelines(filtered)
            written += 1
    return written


def run(cmd, log_path):
    t0 = time.time()
    with open(log_path, 'w') as log:
        r = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
    return r.returncode, time.time() - t0


def parse_eval(stdout):
    m = re.search(r'recall=([\d.]+)%.*prec=([\d.]+)%', stdout)
    if not m:
        return None
    return float(m.group(1)), float(m.group(2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--force', action='store_true')
    args = ap.parse_args()
    WORK.mkdir(parents=True, exist_ok=True)
    rows = []

    for level in LEVELS:
        tag = f'sp{level}'
        blocks = WORK / f'blocks_{tag}.fa'
        mers = WORK / f'mers_{tag}.tsv'
        val = WORK / f'validate_{tag}'
        contigs = WORK / f'contigs_{tag}.fa'
        blk = HERE / 'stats' / f'msa_{tag}.blocks'
        audit = HERE / 'stats' / f'msa_{tag}.audit'
        log = WORK / f'{tag}.log'

        if args.force or not blocks.exists():
            n = subset_blocks(level, blocks)
            print(f'{tag}: subset MSA blocks={n}', flush=True)

        if args.force or not mers.exists():
            t0 = time.time()
            code, _ = run([str(MERTABLE), str(blocks), '-k', '13', '-c', '4',
                           '-o', str(mers)], log)
            if code != 0:
                print(f'{tag}: mertable failed', flush=True)
                continue
            print(f'{tag}: table built in {time.time()-t0:.1f}s '
                  f'({mers.stat().st_size/1e6:.0f} MB)', flush=True)

        if args.force or not contigs.exists():
            if val.exists():
                shutil.rmtree(val)
            code, dt = run([str(VALIDATE),
                            str(BASE / 'reads' / '5x' / 'reads_R1.fq.gz'),
                            str(BASE / 'reads' / '5x' / 'reads_R2.fq.gz'),
                            '--mers', str(mers), '--type', 'fastq', '-t', '8',
                            '--pigz', str(PIGZ),
                            '--rescue-panels', str(blocks),
                            '--rescue-tau', '0.90', '--t3', '1.5,5,25',
                            '--vote-frac', '0.05', '-o', str(val)], log)
            if code != 0:
                print(f'{tag}: validate failed', flush=True)
                continue
            print(f'{tag}: validate {dt:.0f}s', flush=True)

            code, dt = run([str(ASSEMBLE), str(val), '--mers', str(mers),
                            '--trim', '--panels', str(blocks),
                            '--rescue-tau', '0.65', '--rescue-min-span', '30',
                            '--min-rescue-score', '1.00', '-t', '4',
                            '-o', str(contigs)], log)
            if code != 0:
                print(f'{tag}: assemble failed', flush=True)
                continue
            print(f'{tag}: assemble {dt:.0f}s', flush=True)

        r = subprocess.run([str(EVAL), str(REFS), str(contigs), '0.95', '0.30',
                            str(blk), str(audit)],
                           capture_output=True, text=True)
        ev = parse_eval(r.stdout)
        if ev is None:
            print(f'{tag}: eval failed: {r.stderr[:200]}', flush=True)
            continue
        rec, prec = ev
        mers_n = sum(1 for _ in open(mers)) if mers.exists() else 0
        rows.append({
            'species': level,
            'mers': mers_n,
            'table_mb': f'{mers.stat().st_size/1e6:.1f}',
            'recall': f'{rec:.2f}',
            'precision': f'{prec:.2f}',
        })
        print(f'{tag}: recall={rec:.2f} precision={prec:.2f} '
              f'mers={mers_n}', flush=True)

    cols = ['species', 'mers', 'table_mb', 'recall', 'precision']
    with open(OUT, 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for r in rows:
            f.write('\t'.join(str(r[c]) for c in cols) + '\n')
    print(f'Saved: {OUT}')


if __name__ == '__main__':
    main()
