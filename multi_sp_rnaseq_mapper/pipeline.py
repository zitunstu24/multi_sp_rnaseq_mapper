"""Resumable multi-species RNA-seq workflow. Run: python pipeline.py --help"""
import argparse
import json
import logging
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd
import yaml

LOG = logging.getLogger('rnaseq')
STEPS = ['gtf', 'index', 'download', 'star', 'kallisto', 'combine', 'degs', 'cleanup']

def run(cmd):
    LOG.info('RUN %s', ' '.join(map(str, cmd)))
    subprocess.run(list(map(str, cmd)), check=True)

def valid(path):
    p = Path(path)
    return p.is_file() and p.stat().st_size > 0

def needs(outputs, inputs=(), force=False):
    if force or not all(valid(x) for x in outputs):
        return True
    output_time = min(Path(x).stat().st_mtime for x in outputs)
    return any(Path(x).stat().st_mtime > output_time for x in inputs if Path(x).exists())

def atomic_csv(frame, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode='w', dir=path.parent, suffix='.tmp', delete=False) as f:
        tmp = Path(f.name)
    try:
        frame.to_csv(tmp, index=False)
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)

def fastqs(out, sample, layout):
    return [out / f'{sample}_1.fastq', out / f'{sample}_2.fastq'] if layout == 'paired' else [out / f'{sample}.fastq']

def load(cfg):
    with open(cfg) as f:
        config = yaml.safe_load(f)
    root = Path(config['directories']['output_dir']).expanduser().resolve()
    genomes = Path(config['directories']['genome_dir']).expanduser().resolve()
    metadata = Path(config['master_table']).expanduser().resolve()
    df = pd.read_csv(metadata, dtype={'SRA_ID': str, 'species': str, 'layout': str, 'groups': str, 'respective_controls': str}).fillna('')
    required = {'SRA_ID', 'species', 'layout', 'groups', 'respective_controls'}
    if not required.issubset(df.columns):
        raise ValueError(f'Missing master table columns: {required - set(df.columns)}')
    if df.SRA_ID.duplicated().any():
        raise ValueError('Duplicate SRA_ID values in master table')
    if not df.layout.str.lower().isin(['single', 'paired']).all():
        raise ValueError('layout must be single or paired')
    return config, root, genomes, metadata, df

def prepare_gtf(genomes, species, force):
    gtf = genomes / f'{species}.gtf'
    gff = genomes / f'{species}.gff3'
    if gff.exists() and needs([gtf], [gff], force):
        run(['gffread', gff, '-T', '-o', gtf])
    if not valid(gtf):
        raise FileNotFoundError(f'GTF missing: {gtf}; expected {gff} for conversion')
    return gtf

def star_index(genomes, out, species, threads, force):
    fasta, gtf = genomes / f'{species}_genome.fa', genomes / f'{species}.gtf'
    directory = out / f'{species}_genome_index'
    directory.mkdir(parents=True, exist_ok=True)
    marker = directory / '.pipeline_complete.json'
    essentials = [directory / x for x in ['Genome', 'SA', 'SAindex', 'genomeParameters.txt']]
    if needs(essentials + [marker], [fasta, gtf], force):
        marker.unlink(missing_ok=True)
        run(['STAR', '--runThreadN', threads, '--runMode', 'genomeGenerate', '--genomeDir', directory,
             '--genomeFastaFiles', fasta, '--sjdbGTFfile', gtf, '--sjdbOverhang', '99'])
        if not all(valid(x) for x in essentials):
            raise RuntimeError(f'Incomplete STAR index for {species}')
        marker.write_text(json.dumps({'fasta': str(fasta), 'gtf': str(gtf)}))
    else:
        LOG.info('SKIP STAR index %s', species)
    return directory

def kallisto_index(genomes, species, force):
    candidates = [genomes / f'{species}_cds.fa', genomes / f'{species}_transcripts.fa']
    candidates += sorted(genomes.glob(f'{species}*cds.fa'))
    source = next((p for p in candidates if valid(p)), None)
    if source is None:
        raise FileNotFoundError(f'No CDS/transcript FASTA found for {species}')
    index = genomes / f'{species}.idx'
    if needs([index], [source], force):
        index.unlink(missing_ok=True)
        run(['kallisto', 'index', '-i', index, source])
    else:
        LOG.info('SKIP Kallisto index %s', species)
    return index

def download(out, sample, layout, force):
    paths = fastqs(out, sample, layout)
    if not needs(paths, force=force):
        LOG.info('SKIP download %s', sample)
        return
    for p in paths:
        p.unlink(missing_ok=True)
    run(['fasterq-dump', '--split-files', '-O', out, sample] if layout == 'paired' else ['fasterq-dump', '-O', out, sample])
    if not all(valid(p) for p in paths):
        raise RuntimeError(f'Download incomplete: {sample}')

def star(out, species, sample, layout, index, threads, force, stranded):
    directory = out / species / 'star' / sample
    directory.mkdir(parents=True, exist_ok=True)
    result = directory / 'ReadsPerGene.out.tab'
    reads = fastqs(out, sample, layout)
    if not needs([result], reads, force):
        LOG.info('SKIP STAR %s', sample)
        return
    if not all(valid(p) for p in reads):
        raise FileNotFoundError(f'FASTQ missing for STAR {sample}; rerun download or restore reads')
    result.unlink(missing_ok=True)
    run(['STAR', '--runThreadN', threads, '--genomeDir', index, '--readFilesIn', *reads,
         '--outFileNamePrefix', str(directory) + '/', '--outSAMtype', 'None', '--quantMode', 'GeneCounts'])
    if not valid(result):
        raise RuntimeError(f'STAR output missing: {result}')

def kallisto(out, species, sample, layout, index, threads, force, frag_length, frag_sd):
    directory = out / species / 'kallisto' / sample
    directory.mkdir(parents=True, exist_ok=True)
    result = directory / 'abundance.tsv'
    reads = fastqs(out, sample, layout)
    if not needs([result], reads, force):
        LOG.info('SKIP Kallisto %s', sample)
        return
    if not all(valid(p) for p in reads):
        raise FileNotFoundError(f'FASTQ missing for Kallisto {sample}; rerun download or restore reads')
    result.unlink(missing_ok=True)
    cmd = ['kallisto', 'quant', '-i', index, '-o', directory, '-t', threads]
    if layout == 'single':
        cmd += ['--single', '-l', frag_length, '-s', frag_sd]
    run(cmd + reads)
    if not valid(result):
        raise RuntimeError(f'Kallisto output missing: {result}')

def combine(out, species, samples, stranded, force):
    for method in ('star', 'kallisto'):
        target = out / f'combined_{method}_counts_{species}.csv'
        paths = [out / species / method / s / ('ReadsPerGene.out.tab' if method == 'star' else 'abundance.tsv') for s in samples]
        if not all(valid(p) for p in paths):
            raise FileNotFoundError(f'Missing {method} sample outputs: {[str(p) for p in paths if not valid(p)]}')
        if not needs([target], paths, force):
            LOG.info('SKIP combine %s %s', species, method)
            continue
        combined = None
        for s, path in zip(samples, paths):
            if method == 'star':
                data = pd.read_csv(path, sep='\t', header=None).iloc[4:, [0, {'unstranded': 1, 'forward': 2, 'reverse': 3}[stranded]]]
                data.columns = ['gene', s]
            else:
                data = pd.read_csv(path, sep='\t')[['target_id', 'est_counts', 'tpm']]
                data.columns = ['gene', f'{s}_counts', f'{s}_tpm']
            if data.gene.duplicated().any():
                raise ValueError(f'Duplicate identifiers in {path}')
            combined = data if combined is None else combined.merge(data, on='gene', how='outer', validate='one_to_one')
        if combined.isna().any().any():
            raise ValueError(f'Inconsistent feature identifiers among {method} samples for {species}; refusing to fill missing features with zero')
        atomic_csv(combined, target)
        LOG.info('WROTE %s', target)

def degs(out, species, metadata, threads, force, methods, strandedness='unstranded'):
    """Delegate expression matrix creation and DE testing to DEGs.py."""
    from .DEGs import run_degs
    run_degs(out, species, metadata, threads=threads, force=force,
             methods=methods, strandedness=strandedness)

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', required=True)
    parser.add_argument('--from-step', choices=STEPS, default='gtf')
    parser.add_argument('--to-step', choices=STEPS, default='cleanup')
    parser.add_argument('--only-step', choices=STEPS)
    parser.add_argument('--skip', nargs='*', choices=STEPS, default=[])
    parser.add_argument('--force', action='store_true', help='Recompute selected stages')
    parser.add_argument('--species', nargs='+', help='Process only specified species')
    parser.add_argument('--samples', nargs='+', help='Process only specified SRA IDs')
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument('--strandedness', choices=['unstranded', 'forward', 'reverse'], default='unstranded')
    parser.add_argument('--fragment-length', type=float, default=150)
    parser.add_argument('--fragment-sd', type=float, default=20)
    parser.add_argument('--methods', nargs='+', choices=['star', 'kallisto'], default=['star', 'kallisto'])
    parser.add_argument('--cleanup-fastq', action='store_true', help='Explicit opt-in; delete reads only after both quantifications exist')
    parser.add_argument('--dry-run', action='store_true', help='Print stage selection without executing')
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
    if args.threads < 1:
        parser.error('--threads must be >= 1')
    if args.only_step:
        selected = [args.only_step]
    else:
        if STEPS.index(args.from_step) > STEPS.index(args.to_step):
            parser.error('--from-step must precede --to-step')
        selected = STEPS[STEPS.index(args.from_step):STEPS.index(args.to_step)+1]
    selected = [s for s in selected if s not in args.skip]
    if not args.cleanup_fastq:
        selected = [s for s in selected if s != 'cleanup']
    cfg, out, genomes, master, df = load(args.config)
    if args.species:
        df = df[df.species.isin(args.species)]
    if args.samples:
        df = df[df.SRA_ID.isin(args.samples)]
    if df.empty:
        parser.error('No samples selected')
    if args.dry_run:
        print('Selected steps:', ', '.join(selected))
        print('Samples:', ', '.join(df.SRA_ID))
        return
    out.mkdir(parents=True, exist_ok=True)
    for species in df.species.unique():
        sub = df[df.species == species]
        if 'gtf' in selected:
            prepare_gtf(genomes, species, args.force)
        if 'index' in selected:
            if 'star' in args.methods:
                star_index(genomes, out, species, args.threads, args.force)
            if 'kallisto' in args.methods:
                kallisto_index(genomes, species, args.force)
        for row in sub.itertuples(index=False):
            sample, layout = row.SRA_ID, row.layout.lower()
            if 'download' in selected:
                download(out, sample, layout, args.force)
            if 'star' in selected and 'star' in args.methods:
                star(out, species, sample, layout, out / f'{species}_genome_index', args.threads, args.force, args.strandedness)
            if 'kallisto' in selected and 'kallisto' in args.methods:
                kallisto(out, species, sample, layout, genomes / f'{species}.idx', args.threads, args.force, args.fragment_length, args.fragment_sd)
            if 'cleanup' in selected:
                star_file = out / species / 'star' / sample / 'ReadsPerGene.out.tab'
                kallisto_file = out / species / 'kallisto' / sample / 'abundance.tsv'
                if not all(valid(p) for p in (star_file, kallisto_file)):
                    raise RuntimeError(f'Refusing to delete FASTQs: incomplete quantification for {sample}')
                for path in fastqs(out, sample, layout):
                    path.unlink(missing_ok=True)
                LOG.info('Cleaned FASTQs for %s', sample)
        if 'combine' in selected:
            combine(out, species, sub.SRA_ID.tolist(), args.strandedness, args.force)
        if 'degs' in selected:
            degs(out, species, sub, args.threads, args.force, args.methods, args.strandedness)
    LOG.info('Pipeline finished successfully')

if __name__ == '__main__':
    main()
