# Resumable multi-species RNA-seq pipeline (review edition)

Run `python pipeline.py --config path/to/config.yml --dry-run` to inspect stage selection.

Dependencies: Python with pandas and PyYAML; STAR, kallisto, fasterq-dump, gffread in PATH; PyDESeq2 for `degs`.

The YAML configuration must provide `directories.output_dir`, `directories.genome_dir`, and `master_table`. The master CSV must contain `SRA_ID`, `species`, `layout`, `groups`, `respective_controls`. Expected references: `<species>_genome.fa`, `<species>.gtf` or `<species>.gff3`, and `<species>_cds.fa` (or `<species>_transcripts.fa`).

Examples:

```
python pipeline.py --config config/default_config.yml
python pipeline.py --config config/default_config.yml --from-step star --to-step combine
python pipeline.py --config config/default_config.yml --only-step degs
python pipeline.py --config config/default_config.yml --only-step star --force
python pipeline.py --config config/default_config.yml --only-step kallisto --samples SRR123456
python pipeline.py --config config/default_config.yml --skip kallisto --to-step star
```

## Important limitations / validation before publication

- This is a code-level revision and passed Python syntax compilation; no SRA datasets or bioinformatics binaries were available for an end-to-end test.
- Existing outputs are reused when nonempty (and not older than applicable immediate inputs). This is **not** full provenance validation: changing tool versions, options, reference content with unchanged timestamps, or sample metadata may require `--force`.
- STAR index created by the new script uses a marker file. Pre-existing STAR indices without that marker are rebuilt.
- STAR counts use `--strandedness unstranded` by default; set forward/reverse based on experimental library prep, NOT paired-end status.
- Kallisto `est_counts` are rounded for PyDESeq2. This is a pragmatic approximation, not transcript-to-gene summarization. For transcript-level FASTA, aggregate with a validated transcript-to-gene mapping (e.g. tximport workflow) before gene-level DEGs.
- Single-end Kallisto requires accurate fragment length and SD (`--fragment-length`, `--fragment-sd`); defaults of 150/20 are inherited assumptions, not validated sample-specific measurements.
- STAR gene-counts and Kallisto estimates should not be compared as identical quantification methods.
- DESeq2 results are written separately to `output/degs/<species>/<method>/<treatment>_vs_<control>.csv`; the original combined DEG-wide CSV is not generated.
- `--samples` + `--only-step combine` will rebuild matrices from only the selected samples; do not use a sample filter when creating final species-wide matrices.
- The pipeline will not automatically redownload FASTQs if a downstream-only run needs files that were cleaned. Include the download stage, or restore the FASTQs.
- Cleanup is opt-in and per sample, and requires both STAR and Kallisto output files. It does not verify data integrity beyond nonempty files.
- Species, treatment and control strings are used as path segments: use safe names without `/` or `..`.
- For stronger provenance and HPC scheduling, a Snakemake/Nextflow DAG with checksums and environment locking is recommended.
