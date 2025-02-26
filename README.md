
# Cross Species SRA Download, Alignment and Benchmarking Tool

This tool automates the process of downloading SRA samples, aligning them to a reference genome using the STAR aligner, and cleaning up fastq files afterward. The tool is designed to be modular, configurable, and easy to install, making it suitable for bioinformatics workflows.

## Features

Download SRA Samples: 

Automatically download SRA samples using the fasterq-dump tool from the SRA Toolkit.

Genome Indexing: 
Index genome FASTA files with STAR for efficient alignments.

Run STAR Alignment: 

Align downloaded fastq files to species-specific genomes using STAR.

Cleanup: 
Deletes fastq files after alignment to save space.

Configurable: 
Use a YAML configuration file to define directories, output locations, and other settings.

## Installation

Ensure that the following tools are installed on your system: 

SRA Toolkit: Used to download SRA files.
Installation instructions: SRA Toolkit
STAR: Required for genome indexing and read alignment.
Installation instructions: STAR Aligner
Python 3.8+: The tool requires Python 3.8 or higher.

## How to use it

```bash
git clone Git link
cd multi_sp_rnaseq_mapper
pip install .
```

Configuration
This tool uses a YAML file for configuration. You can customize various aspects, such as directories for genomes, output locations, and STAR parameters.

An example configuration file can be found at `config/default_config.yml`. Here's a typical structure:


```yml
# config/default_config.yml

# Configuration 
directories:
  genome_dir: "/winmounts/khayer6872/AG-Gutjahr/Group_Projects/abul_k/p0424-AK/genome_gff3"
  output_dir: "./data/outputs/"


master_table: "./data/inputs/master_table.csv" 

star:
  threads: 4

```

### Important Directories

`genome_dir`: Path to the directory containing the species' genome FASTA and GFF3 files.
`output_dir`: Path where outputs (such as BAM files, read_counts etc) will be stored.

Usage
After installation, you can use the tool by providing the path to the configuration file.

### Master Table
The master table should be a CSV file with the following format:

```csv
SRA_ID,species,layout
SRR123456,rice,paired
SRR789101,brachy,single
SRR112233,lotus,paired

```

## Running the Tool

You can run the SRA Alignment Tool using the following command:
```bash
RnaMapper --config config/default_config.yml

```
### Command-line Arguments

```bash

--config: Path to the configuration YAML file.

```

## This will:

1. Download each SRA sample listed in master_table.csv.
2. Align the fastq files to the respective species' genome using STAR.
3. Delete the fastq files once alignment is complete.
4. Align all samples per species in the output directory
5. Analyse differentilly expressed genes with designed samples

## Contributing
Feel free to submit issues or pull requests to improve the functionality of the tool. Contributions are welcome!

## Contact
Feel free to contact Abul KHayer(MPI-MP) : abul.khayer@mpimp-golm.mpg.de or zitunstu24@gmail.com
