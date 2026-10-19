# Multi-species RNA-seq pipeline — restartable version

This project provides a Python CLI for acquiring SRA reads, preparing reference indices, counting reads with STAR, estimating abundance with Kallisto, combining per-species count tables, and calling differential expression with PyDESeq2. It is intended for species with separate reference genomes and annotations.

**Status:** Syntax-tested and packaging-tested; no end-to-end RNA-seq execution has been performed. STAR, Kallisto, SRA Toolkit and gffread are external executables and **are not installed by pip**.

## 1. Folder layout

```text
project/
├── pyproject.toml
├── setup.py
├── pipeline.py                   # source-checkout convenience runner
├── multi_sp_rnaseq_mapper/
│   ├── __init__.py
│   ├── __main__.py
│   ├── pipeline.py               # actively maintained restartable CLI
│   └── ...                       # original helper modules, retained
├── config/
│   └── default_config.yml       # create from example below
├── metadata/
│   └── samples.csv              # create from example below
├── references/
│   ├── rice_genome.fa
│   ├── rice.gff3                 # or rice.gtf
│   └── rice_cds.fa               # or rice_transcripts.fa
└── results/                      # created on execution
```

## 2. Installation

Python **3.10+** is required. From the project root (where `pyproject.toml` lives):

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e '.[degs]'
```

If you do not need the DEGs stage, `python -m pip install -e .` installs only pandas and PyYAML. External programs must be installed separately in your HPC environment or Conda environment:

```bash
command -v STAR
command -v kallisto
command -v fasterq-dump
command -v gffread
```

Check versions with `STAR --version`, `kallisto version`, `fasterq-dump --version`, and `gffread --version`. Use a tested environment lock file for reproducibility. PyDESeq2 imports are deferred until running `degs`.

After installation, these commands invoke **the same pipeline**:

```bash
RnaMapper --help
multi-sp-rnaseq-mapper --help
python -m multi_sp_rnaseq_mapper --help
python pipeline.py --help          # in the source checkout only
python main.py --help              # legacy alias in the source checkout
```

`setup.py` is a compatibility shim: dependency metadata, version and console scripts are defined in `pyproject.toml`. Do not maintain two conflicting dependency lists.

## 3. Input configuration

Create `config/default_config.yml`:

```yaml
directories:
  genome_dir: /absolute/path/to/project/references
  output_dir: /absolute/path/to/project/results
master_table: /absolute/path/to/project/metadata/samples.csv
```

The current implementation resolves YAML path values relative to the **current working directory**, not relative to the YAML file. Absolute paths are recommended.

Create the master CSV with the following **exact headers**:

```csv
SRA_ID,species,layout,groups,respective_controls
SRR000001,rice,paired,control,
SRR000002,rice,paired,control,
SRR000003,rice,paired,treatment,control
SRR000004,rice,paired,treatment,control
```

The SRA identifiers above are placeholders, **not asserted to be valid accessions**. Each `SRA_ID` must be unique. `species` must match reference filenames; `layout` is `single` or `paired`; `groups` is the experimental group; `respective_controls` on treatment samples gives the exact corresponding control group name. Control samples generally leave that field empty. Comparisons require at least two biological replicates in each group according to the current checks; greater replication may be advisable for real analyses. Ensure group names do not contain slashes or other unsafe path characters.

Each species should have:

- `<species>_genome.fa` (genomic FASTA for STAR).
- `<species>.gtf`, or `<species>.gff3` convertible via gffread.
- `<species>_cds.fa` or `<species>_transcripts.fa` as the Kallisto index reference. The implementation additionally searches `<species>*cds.fa`.

**Important:** CDS/transcript features in Kallisto abundance are not automatically aggregated into gene counts. Do not interpret Kallisto DEGs as gene-level results unless your FASTA target identifiers truly represent genes, or you have implemented and validated transcript-to-gene aggregation.

## 4. Pipeline steps and outputs

Stages execute in this fixed order:

| Stage | Task | Principal output / reuse check |
|---|---|---|
| `gtf` | Convert GFF3 to GTF if necessary | `references/<species>.gtf` |
| `index` | Index STAR genome and Kallisto target FASTA | `results/<species>_genome_index/` and `references/<species>.idx` |
| `download` | Download SRA FASTQs | `results/<SRA_ID>.fastq` or `_1.fastq` and `_2.fastq` |
| `star` | STAR alignment with gene-count output | `results/<species>/star/<SRA_ID>/ReadsPerGene.out.tab` |
| `kallisto` | Quantification with Kallisto | `results/<species>/kallisto/<SRA_ID>/abundance.tsv` |
| `combine` | Join counts by gene/target ID per species | `results/combined_star_counts_<species>.csv`, `results/combined_kallisto_counts_<species>.csv` |
| `degs` | Differential expression for each treatment vs control | `results/degs/<species>/<method>/<treatment>_vs_<control>.csv` |
| `cleanup` | Optional removal of sample FASTQs | No count outputs removed |

STAR writes gene counts with `--quantMode GeneCounts` and no BAM (`--outSAMtype None`). `--strandedness` selects which of the STAR gene-count columns is used during `combine` (it does not change the STAR counting command). Choose library strandedness using library preparation metadata, **not** the single/paired layout. Kallisto single-end mode uses specified fragment-length distribution parameters.

## 5. Commands

Full run (default skips nonempty, non-stale immediate outputs):

```bash
RnaMapper --config config/default_config.yml
```

Preview selection without invoking external tools:

```bash
RnaMapper --config config/default_config.yml --dry-run
```

Restart from STAR through DEG analysis:

```bash
RnaMapper --config config/default_config.yml --from-step star --to-step degs
```

Run only one stage:

```bash
RnaMapper --config config/default_config.yml --only-step degs
RnaMapper --config config/default_config.yml --only-step combine
```

Limit to a species or SRA sample:

```bash
RnaMapper --config config/default_config.yml --species rice --from-step download --to-step kallisto
RnaMapper --config config/default_config.yml --samples SRR000003 --only-step star
```

Recreate STAR outputs, with an explicitly chosen strandedness for downstream matrix building:

```bash
RnaMapper --config config/default_config.yml --only-step star --force --threads 16
RnaMapper --config config/default_config.yml --only-step combine --strandedness reverse --force
```

Skip stages or choose a method:

```bash
RnaMapper --config config/default_config.yml --skip download --from-step star --to-step kallisto
RnaMapper --config config/default_config.yml --methods star --to-step star
RnaMapper --config config/default_config.yml --methods kallisto --only-step degs
```

Cleanup, **opt-in only**, and requires both quantifications already present:

```bash
RnaMapper --config config/default_config.yml --only-step cleanup --cleanup-fastq
```

After FASTQs are cleaned, a forced re-alignment or re-quantification needs re-downloaded/restored reads. A default full run can also download again because the FASTQs are absent.

## 6. CLI reference

| Argument | Effect |
|---|---|
| `--config PATH` | Required YAML configuration |
| `--from-step STAGE` | First included stage (default `gtf`) |
| `--to-step STAGE` | Last included stage (default `cleanup`, though cleanup disabled unless opted in) |
| `--only-step STAGE` | Run just this stage; supersedes range selection |
| `--skip [STAGE ...]` | Omit listed stages |
| `--force` | Rebuild selected outputs even if present |
| `--species NAME [NAME ...]` | Filter metadata by species |
| `--samples ID [ID ...]` | Filter metadata by SRA IDs |
| `--threads N` | Threads for STAR/Kallisto and DEGs (default 8) |
| `--methods {star,kallisto} ...` | Enable quantification methods (default both) |
| `--strandedness {unstranded,forward,reverse}` | STAR count column for combine (default unstranded) |
| `--fragment-length FLOAT` | Kallisto single-end mean fragment length (default 150) |
| `--fragment-sd FLOAT` | Kallisto single-end fragment length SD (default 20) |
| `--cleanup-fastq` | Allow sample-specific FASTQ deletion at cleanup step |
| `--dry-run` | Print selected stages and samples; **does not** validate all dependencies |

Valid STAGE choices: `gtf`, `index`, `download`, `star`, `kallisto`, `combine`, `degs`, `cleanup`.

## 7. Restarting safely

The pipeline considers a nonempty output reusable if it is not older than any **existing immediate input** checked for that stage. For STAR and Kallisto, when raw FASTQs have been removed, an existing nonempty count output may still be reused. This is practical restart behavior, **not** full checksum/version provenance, and it does not prove results are biologically valid. Use `--force` after reference, parameter or software changes and consider retaining read-only raw data.

`--from-step star` assumes the STAR index exists and will not automatically create it; include `index` when needed. `--only-step degs` assumes the combined CSV exists. `--only-step combine` currently requires **both** STAR and Kallisto sample outputs, even if `--methods` selects only one. The `--methods` flag filters indexing, quantification and DEG steps but not the combine stage. Do not request `combine` in a method-specific run when the other method's outputs are absent.

**Avoid `--samples` with `combine` or `degs`** for final species-wide outputs: a filtered sample list can produce an incomplete matrix, and can overwrite an existing combined CSV if `--force` is used. For a single-sample repair, run only download/STAR/Kallisto for that sample, then rerun `combine` and `degs` for the complete species without `--samples`.

## 8. Differential expression details

For each distinct treatment group with an assigned control, counts are assembled in sample × feature orientation and DESeq2-style statistics are calculated using PyDESeq2 with an explicit treatment-versus-control contrast. The code rounds Kallisto estimated counts to integers, which is an approximation and **not** a substitute for established transcript-level quantification import methods. It does not explicitly adjust for sequencing batch or other covariates. Inspect metadata quality, feature identifiers, replication, count distributions, experimental design, and the output statistics before biological interpretation.

## 9. Troubleshooting

| Symptom | Check |
|---|---|
| `command not found: STAR` / `kallisto` / `fasterq-dump` / `gffread` | Activate the environment and confirm binaries are in `PATH` |
| `ModuleNotFoundError: pydeseq2` | Install `python -m pip install -e '.[degs]'` or install PyDESeq2 in active environment |
| `GTF missing` | Match `species` with `<species>.gff3` or `<species>.gtf` under genome directory |
| `No CDS/transcript FASTA found` | Add matching `species_cds.fa` / `species_transcripts.fa` |
| `FASTQ missing for STAR/Kallisto` | Include `download` or restore reads; avoid cleanup before reprocessing |
| `Missing ... sample outputs` in combine | Run missing STAR/Kallisto sample stages; combine expects both methods |
| `Inconsistent feature identifiers` | Confirm all samples share the same annotation/index and feature ID space |
| `Missing count columns` / control missing | Check master CSV identifiers and control names |
| Index or results unexpectedly reused | Check timestamps; use `--force` after input/parameter changes |
| Installation succeeds but running fails | Python packaging cannot install external genomics binaries; verify tool dependencies separately |

## 10. Reproducibility checklist

Record reference genome/annotation checksums, tool versions, Python environment versions, samples and experimental design, command line and runtime logs; keep metadata stable; verify library strandedness and sequencing layout; check output QC; rerun with `--force` when dependencies change. To obtain complete dependency-aware execution and HPC job scheduling, migrate stage definitions to a workflow engine such as Snakemake or Nextflow.

## 11. Packaging / development

Editable install: `python -m pip install -e '.[degs]'`. Standard install: `python -m pip install .`. Build distributions after installing `build`: `python -m build`. Console script mapping is declared under `[project.scripts]` in `pyproject.toml`. `setup.py` exists for older build tooling and delegates metadata to the same file. The original helper modules are retained for compatibility, but the new CLI runs `multi_sp_rnaseq_mapper/pipeline.py`, not the legacy `process_samples.py` + `DEGs.py` orchestration.
