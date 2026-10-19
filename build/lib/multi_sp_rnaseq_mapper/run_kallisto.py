
import os
import subprocess

def create_kallisto_index(genome_dir):
    """Create Kallisto index from a transcriptome FASTA file."""
    
    # Loop through each cds fasta file in the directory
    for file_name in os.listdir(genome_dir):
        if file_name.endswith("cds.fa"):
            cds_file = os.path.join(genome_dir, file_name)
            index = file_name.split('_')[0]
            index = f"{index}.idx"
            index_file = os.path.join(genome_dir, index)
            if not os.path.exists(index_file):
                try:
                    print("Creating Kallisto index...")
                    command = ['kallisto', 'index', '-i', index_file, cds_file]
                    subprocess.run(command, check=True)
                except subprocess.CalledProcessError as e:
                    print(f"Error creating Kallisto index: {e}")

def run_kallisto(sra_id, species, index, output_dir, layout):
    """Run Kallisto to quantify transcript abundance based on layout (single or paired)."""

    # Paths for fastq files
    fastq1 = os.path.join(output_dir, f"{sra_id}_1.fastq")
    fastq2 = os.path.join(output_dir, f"{sra_id}_2.fastq") if layout == 'paired' else None
    single_fastq = os.path.join(output_dir, f"{sra_id}.fastq") if layout == 'single' else None

    # Directory to save Kallisto output for each species and sample
    species_output_dir = os.path.join(output_dir, species, "kallisto")
    kallisto_output_dir = os.path.join(species_output_dir, sra_id)
    os.makedirs(kallisto_output_dir, exist_ok=True)

    try:
        print(f"Running Kallisto for {sra_id} on species {species} ({layout}-end)...")
        command = ['kallisto', 'quant', '-i', index, '-o', kallisto_output_dir]

        if layout == 'paired':
            command += ['-t', '50', fastq1, fastq2]
        elif layout == 'single':
            command += ['-t', '50', '--single', '-l', '150', '-s', '20', single_fastq]
        else:
            raise ValueError(f"Unknown layout: {layout}")

        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running Kallisto for {sra_id}: {e}")
