"""DEG analysis and expression matrices from SRA-level STAR/Kallisto outputs.

Inputs are metadata with SRA_ID/species/groups/respective_controls and either
per-sample quantification files or previously combined counts matrices.
"""
from __future__ import annotations

import logging
from pathlib import Path
import re

import numpy as np
import pandas as pd

LOG = logging.getLogger('rnaseq')


def _atomic(df: pd.DataFrame, path: Path, *, index: bool = False) -> None:
    import filecmp
    import os
    import tempfile
    if index:
        df = df.reset_index()
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile('w', encoding='utf-8', suffix='.csv.tmp',
                                     dir=path.parent, delete=False) as tmp:
        tmp_path = Path(tmp.name)
        df.to_csv(tmp, index=False)
    try:
        if not (path.exists() and filecmp.cmp(tmp_path, path, shallow=False)):
            os.replace(tmp_path, path)
    finally:
        tmp_path.unlink(missing_ok=True)


def _read_sample(out: Path, species: str, method: str, sra: str, strandedness: str) -> pd.DataFrame:
    root = out / species / method / sra
    if method == 'star':
        path = root / 'ReadsPerGene.out.tab'
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f'STAR counts absent for {sra}: {path}')
        data = pd.read_csv(path, sep='\t', header=None)
        if data.shape[1] < 4:
            raise ValueError(f'Malformed STAR counts: {path}')
        col = {'unstranded': 1, 'forward': 2, 'reverse': 3}[strandedness]
        data = data.iloc[4:, [0, col]].copy()
        data.columns = ['gene', sra]
    elif method == 'kallisto':
        path = root / 'abundance.tsv'
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f'Kallisto abundance absent for {sra}: {path}')
        data = pd.read_csv(path, sep='\t')
        missing = {'target_id', 'est_counts', 'tpm'} - set(data.columns)
        if missing:
            raise ValueError(f'Malformed Kallisto file {path}; missing {missing}')
        data = data[['target_id', 'est_counts', 'tpm']].copy()
        data.columns = ['gene', sra, f'{sra}_tpm']
    else:
        raise ValueError(f'Unknown method: {method}')
    if data['gene'].isna().any() or data['gene'].astype(str).duplicated().any():
        raise ValueError(f'Missing/duplicate gene or target identifiers in {path}')
    data['gene'] = data['gene'].astype(str)
    for col in data.columns[1:]:
        data[col] = pd.to_numeric(data[col], errors='raise')
        if not np.isfinite(data[col].to_numpy(dtype=float)).all() or (data[col] < 0).any():
            raise ValueError(f'Invalid expression values in {path}, column {col}')
    return data.set_index('gene')


def _existing_combined(out: Path, species: str, method: str, sample_ids: list[str]) -> tuple[pd.DataFrame, pd.DataFrame | None]:
    path = out / f'combined_{method}_counts_{species}.csv'
    if not path.is_file():
        raise FileNotFoundError(f'Neither complete per-sample results nor a combined matrix found: {path}')
    df = pd.read_csv(path)
    if 'gene' not in df.columns or df['gene'].isna().any() or df['gene'].duplicated().any():
        raise ValueError(f'Invalid gene/target column in {path}')
    df = df.set_index('gene')
    cols = {s: (s if s in df.columns else f'{s}_counts') for s in sample_ids}
    missing = [s for s, c in cols.items() if c not in df.columns]
    if missing:
        raise ValueError(f'SRA IDs absent from {path}: {missing}')
    raw = df[[cols[s] for s in sample_ids]].copy()
    raw.columns = sample_ids
    tpmcols = [f'{s}_tpm' for s in sample_ids]
    tpm = df[tpmcols].copy() if method == 'kallisto' and all(c in df.columns for c in tpmcols) else None
    if tpm is not None:
        tpm.columns = sample_ids
    return raw, tpm


def _matrix(out: Path, species: str, method: str, ids: list[str], strandedness: str) -> tuple[pd.DataFrame, pd.DataFrame | None]:
    filenames = ['ReadsPerGene.out.tab' if method == 'star' else 'abundance.tsv']
    paths = [out / species / method / s / filenames[0] for s in ids]
    if not all(p.is_file() and p.stat().st_size > 0 for p in paths):
        LOG.info('Some %s sample-level outputs missing; using existing combined matrix for %s', method, species)
        return _existing_combined(out, species, method, ids)
    chunks = [_read_sample(out, species, method, s, strandedness) for s in ids]
    reference = chunks[0].index
    for s, chunk in zip(ids[1:], chunks[1:]):
        if set(chunk.index) != set(reference):
            raise ValueError(f'{method}: gene/target identifiers differ for sample {s}; cannot safely combine')
    full = pd.concat([c.reindex(reference) for c in chunks], axis=1)
    counts = full[ids].copy()
    tpm = full[[f'{s}_tpm' for s in ids]].copy() if method == 'kallisto' else None
    if tpm is not None:
        tpm.columns = ids
    return counts, tpm


def _numeric(df: pd.DataFrame, method: str) -> pd.DataFrame:
    result = df.apply(pd.to_numeric, errors='raise')
    vals = result.to_numpy(dtype=float)
    if not np.isfinite(vals).all() or (vals < 0).any():
        raise ValueError(f'{method}: counts include missing, infinite, or negative values')
    if (result.sum(axis=0) <= 0).any():
        raise ValueError(f'{method}: at least one sample has zero total counts')
    return result


def _scaled(counts: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """log2(CPM+1), followed by row-wise z-scoring (ddof=0)."""
    cpm = counts.div(counts.sum(axis=0), axis=1) * 1_000_000.0
    log2cpm = np.log2(cpm + 1)
    means = log2cpm.mean(axis=1)
    std = log2cpm.std(axis=1, ddof=0)
    z = log2cpm.sub(means, axis=0).div(std.replace(0, np.nan), axis=0).fillna(0.0)
    return log2cpm, z


def _safe_name(label: str) -> str:
    return re.sub(r'[^A-Za-z0-9_.-]+', '_', str(label)).strip('._') or 'unnamed'


def _contrasts(samples: pd.DataFrame, species: str) -> list[tuple[str, str]]:
    contrasts = []
    for treatment in samples['groups'].drop_duplicates():
        vals = samples.loc[samples.groups.eq(treatment), 'respective_controls'].dropna().astype(str).str.strip()
        controls = vals[vals.ne('')].unique().tolist()
        if not controls:
            continue
        if len(controls) != 1:
            raise ValueError(f'{species}/{treatment}: conflicting controls {controls}')
        control = controls[0]
        if control == treatment:
            raise ValueError(f'{species}: treatment and control are identical: {treatment}')
        if control not in set(samples.groups):
            raise ValueError(f'{species}/{treatment}: control group {control!r} absent from metadata')
        contrasts.append((treatment, control))
    return contrasts


def run_degs(out: Path, species: str, metadata: pd.DataFrame, threads: int = 8,
             force: bool = False, methods: list[str] | tuple[str, ...] = ('star', 'kallisto'),
             strandedness: str = 'unstranded') -> None:
    """Write all-sample matrices + per-contrast and combined DEG tables.

    All matrices include all SRA_IDs for the selected species, even controls.
    Z-scores are computed gene-wise across *all* samples of that species, using
    log2(CPM+1). Constant rows receive z=0.
    """
    required = {'SRA_ID', 'species', 'groups', 'respective_controls'}
    if not required.issubset(metadata.columns):
        raise ValueError(f'Metadata is missing columns: {sorted(required - set(metadata.columns))}')
    sub = metadata.loc[metadata.species.astype(str).eq(species), list(required)].copy()
    for col in ('SRA_ID', 'groups', 'respective_controls'):
        sub[col] = sub[col].fillna('').astype(str).str.strip()
    if sub.empty or sub.SRA_ID.eq('').any() or sub.groups.eq('').any() or sub.SRA_ID.duplicated().any():
        raise ValueError(f'{species}: invalid or duplicate SRA_ID / empty groups in metadata')
    ids = sub.SRA_ID.tolist()
    sub = sub.set_index('SRA_ID').loc[ids]
    comparisons = _contrasts(sub, species)
    if not comparisons:
        LOG.warning('%s: no respective_controls contrasts; writing expression matrices only', species)
    dest = out / 'degs' / species
    dest.mkdir(parents=True, exist_ok=True)
    annotation = sub[['groups', 'respective_controls']].reset_index()
    _atomic(annotation, dest / 'sample_metadata.csv')

    for method in methods:
        if method not in ('star', 'kallisto'):
            raise ValueError(f'Unsupported method {method}')
        method_dest = dest / method
        method_dest.mkdir(parents=True, exist_ok=True)
        counts, tpm = _matrix(out, species, method, ids, strandedness)
        counts = _numeric(counts.loc[:, ids], method)
        counts.index.name = 'gene'
        # Store expression and annotations separately: gene rows x SRA_ID columns.
        _atomic(counts, method_dest / 'all_samples_counts.csv', index=True)
        # Compatibility file for downstream scripts (SRA-named columns).
        _atomic(counts, out / f'combined_{method}_counts_{species}.csv', index=True)
        if tpm is not None:
            _atomic(tpm, method_dest / 'all_samples_tpm.csv', index=True)
        log2cpm, z = _scaled(counts)
        _atomic(log2cpm, method_dest / 'all_samples_log2cpm.csv', index=True)
        _atomic(z, method_dest / 'all_samples_zscores.csv', index=True)
        LOG.info('%s/%s: wrote all-sample matrices: %d genes/targets x %d SRA samples', species, method, len(counts), len(ids))

        combined = None
        for treatment, control in comparisons:
            group_sub = sub.loc[sub.groups.isin((treatment, control)), ['groups']].copy()
            found = set(group_sub.groups)
            if found != {treatment, control}:
                raise ValueError(f'{species}: expected {treatment!r}, {control!r}, found {sorted(found)}')
            if (group_sub.groups.value_counts() < 2).any():
                raise ValueError(f'{species}: fewer than 2 replicates in {treatment} or {control}')
            selected_ids = group_sub.index.tolist()
            # Counts are SRA-aligned exactly to metadata index, NEVER aligned by row order.
            x = counts.loc[:, selected_ids].T
            # Round only for DESeq2 model; retain floating-point Kallisto estimates in CSVs.
            x = np.floor(x + 0.5).astype('int64')
            x = x.loc[:, x.sum(axis=0).gt(0)]
            if x.empty:
                raise ValueError(f'{species}/{method}: all features have zero counts')
            group_sub['groups'] = pd.Categorical(group_sub.groups.astype(str), categories=[control, treatment])
            name = f'{_safe_name(treatment)}_vs_{_safe_name(control)}'
            target = method_dest / f'{name}.csv'
            from .pipeline import needs
            source = out / f'combined_{method}_counts_{species}.csv'
            if needs([target], [source], force=force):
                from pydeseq2.dds import DeseqDataSet
                from pydeseq2.ds import DeseqStats
                from pydeseq2.default_inference import DefaultInference
                inference = DefaultInference(n_cpus=threads)
                dds = DeseqDataSet(counts=x, metadata=group_sub, design_factors='groups',
                                   refit_cooks=True, inference=inference)
                dds.deseq2()
                # PyDESeq2 may sanitize group levels (e.g. NAC_OE -> NAC-OE).
                # Derive the fitted levels from the exact SRA-matched observations,
                # not from an assumed string replacement or sample row order.
                observed = dds.obs['groups'].astype(str)
                expected = group_sub['groups'].astype(str)
                if not observed.index.is_unique or set(observed.index) != set(expected.index):
                    raise ValueError(
                        f'{species}/{method}: PyDESeq2 observation IDs no longer match SRA metadata. '
                        f'Expected={expected.index.tolist()}, observed={observed.index.tolist()}')
                label_map = {}
                for original in (control, treatment):
                    sample_ids = expected.index[expected == original]
                    actual_levels = set(observed.loc[sample_ids])
                    if len(actual_levels) != 1:
                        raise ValueError(
                            f'{species}/{method}: group {original!r} maps to multiple PyDESeq2 '
                            f'levels {sorted(actual_levels)}; sample alignment must be checked')
                    label_map[original] = next(iter(actual_levels))
                if label_map[control] == label_map[treatment]:
                    raise ValueError(
                        f'{species}/{method}: PyDESeq2 group-name sanitization collapsed '
                        f'{control!r} and {treatment!r} to the same level')
                LOG.info('DESeq2 contrast: %s (%s) vs %s (%s)',
                         treatment, label_map[treatment], control, label_map[control])
                stats = DeseqStats(dds, contrast=['groups', label_map[treatment],
                                                 label_map[control]], inference=inference)
                stats.summary()
                results = stats.results_df.copy()
                results.index.name = 'gene'
                _atomic(results, target, index=True)
                LOG.info('WROTE DEGs %s: %s vs %s', target, treatment, control)
            else:
                LOG.info('SKIP existing DEGs %s', target)
            result = pd.read_csv(target)
            if 'gene' not in result.columns:
                raise ValueError(f'Missing gene column in DEG result {target}')
            suffix = f'__{name}'
            result = result.rename(columns={c: f'{c}{suffix}' for c in result.columns if c != 'gene'})
            if combined is None:
                combined = result
            else:
                combined = combined.merge(result, on='gene', how='outer', validate='one_to_one')
        if combined is not None:
            _atomic(combined, method_dest / 'all_contrasts_degs.csv')
            _atomic(combined, out / f'combined_DEGs_{method}_{species}.csv')


def get_DEGs(config_file: str, threads: int = 8, force: bool = False,
             methods: tuple[str, ...] = ('star', 'kallisto'), strandedness: str = 'unstranded') -> None:
    """Backward-compatible entry point used by older code."""
    from .pipeline import load
    _, out, _, _, metadata = load(config_file)
    for species in metadata.species.unique():
        run_degs(out, species, metadata, threads=threads, force=force,
                 methods=methods, strandedness=strandedness)
