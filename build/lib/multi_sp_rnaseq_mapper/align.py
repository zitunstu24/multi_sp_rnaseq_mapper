
import os
import subprocess

def run_star(sra_id, species, genome_dir, output_dir, layout, gtf_file=None):
    """Run STAR to align fastq files to the genome based on layout (single or paired)."""
    
    # Prepare paths for fastq files
    fastq1 = os.path.join(output_dir, f"{sra_id}_1.fastq")
    fastq2 = os.path.join(output_dir, f"{sra_id}_2.fastq") if layout == 'paired' else None
    single_fastq = os.path.join(output_dir, f"{sra_id}.fastq") if layout == 'single' else None

    # Create species-specific directory in output folder
    species_output_dir = os.path.join(output_dir, species, "star")
    os.makedirs(species_output_dir, exist_ok=True)

    # STAR output will go into species subdirectory
    star_output_dir = os.path.join(species_output_dir, sra_id)
    os.makedirs(star_output_dir, exist_ok=True)


    try:
        print(f"Running STAR for {sra_id} on species {species} ({layout}-end)...")

        command = ['STAR', '--runThreadN', '100', '--genomeDir', genome_dir, '--outFileNamePrefix', os.path.join(star_output_dir, '')]

        # Adjust STAR command for single or paired layout
        if layout == 'paired':
            command += ['--readFilesIn', fastq1, fastq2]
        elif layout == 'single':
            command += ['--readFilesIn', single_fastq]
        else:
            raise ValueError(f"Unknown layout: {layout}")

        command += ['--outSAMtype', 'None']
        command += ['--quantMode', 'GeneCounts']

        # Include GFF3 file if provided
        if gtf_file:
            command += ['--sjdbGTFfile', gtf_file]

        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running STAR for {sra_id}: {e}")
