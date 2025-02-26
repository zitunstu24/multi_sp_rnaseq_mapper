
import os
import subprocess

def index_genome(species, genome_file, gtf_file, output_dir):
    """Create a STAR genome index for the species if not already exists."""
    index_dir = os.path.join(output_dir, f"{species}_genome_index")
    os.makedirs(index_dir, exist_ok=True)

    # Check if the index already exists
    if not os.listdir(index_dir):  # If directory is empty
        print(f"Indexing genome for {species}...")
        try:
            command = [
                'STAR', '--runThreadN', '100',
                '--runMode', 'genomeGenerate',
                '--genomeDir', index_dir,
                '--genomeFastaFiles', genome_file,
                '--sjdbOverhang', '99',
                '--genomeSAindexNbases', '13'
            ]
            
            # Include the GFF3 file if present
            if gtf_file:
                command += ['--sjdbGTFfile', gtf_file]

            subprocess.run(command, check=True)
            print(f"Genome index created for {species} in {index_dir}.")
        except subprocess.CalledProcessError as e:
            print(f"Error indexing genome for {species}: {e}")
    else:
        print(f"Genome index for {species} already exists at {index_dir}.")
    
    return index_dir
