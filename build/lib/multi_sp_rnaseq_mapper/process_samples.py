from multi_sp_rnaseq_mapper.config import Config
from multi_sp_rnaseq_mapper.download import download_sra
from multi_sp_rnaseq_mapper.align import run_star
from multi_sp_rnaseq_mapper.run_kallisto import create_kallisto_index, run_kallisto
from multi_sp_rnaseq_mapper.genome_index import index_genome
from multi_sp_rnaseq_mapper.cleanup import delete_fastq
from multi_sp_rnaseq_mapper.gff2gtf import convert_gff3_to_gtf
import pandas as pd
import os

def process_samples(config_file):
    # Load the configuration
    config = Config(config_file)

    genome_dir = config.get('directories', 'genome_dir')
    output_dir = config.get('directories', 'output_dir')
    master_table = config.get_master_table()
    gff3_dir = genome_dir

    # Convert GFF3 to GTF format if needed
    convert_gff3_to_gtf(gff3_dir)

    # Prepare Kallisto index if not present

    create_kallisto_index(genome_dir)

    # Read the master table
    df = pd.read_csv(master_table)

    # Loop through each sample in the master table
    for index, row in df.iterrows():
        sra_id = row['SRA_ID']
        species = row['species']
        layout = row['layout'].lower()  # single or paired

        # Define genome and GFF3 file paths
        species_genome = os.path.join(genome_dir, f"{species}_genome.fa")
        species_gtf = os.path.join(genome_dir, f"{species}.gtf")

        # Index the genome if necessary
        index_dir = index_genome(species, species_genome, species_gtf, output_dir)
        
        # Download the SRA sample based on the layout
        download_sra(sra_id, output_dir, layout)
        
        # Run STAR alignment
        run_star(sra_id, species, index_dir, output_dir, layout, species_gtf)
        
        # Run Kallisto quantification
        index = f"{species}.idx"
        index_file = os.path.join(genome_dir, index)
        run_kallisto(sra_id, species, index_file, output_dir, layout)

        # Cleanup fastq files
        delete_fastq(output_dir)

    # Post-processing: Combine counts for STAR and Kallisto as needed
    for species in df['species'].unique():
        print(f"Started combining for {species}")
        species_output_dir = os.path.join(output_dir, species)
        
        # Separate dataframes for STAR and Kallisto
        star_dataframes = pd.DataFrame()
        kallisto_dataframes = pd.DataFrame()
        
        # Process STAR output
        star_output_dir = os.path.join(species_output_dir, "star")
        for subdir, _, files in os.walk(star_output_dir):
            for file in files:
                if file.endswith("ReadsPerGene.out.tab"):
                    file_path = os.path.join(subdir, file)
                    counts_table = pd.read_csv(file_path, header=None, delimiter="\t")
                    
                    # Extract SRA ID from the folder name
                    SRR = os.path.basename(subdir)
                    
                    # Get layout type for determining which column to select
                    layout_check = df[df['SRA_ID'] == SRR]['layout'].iloc[0]
                    
                    # Select columns based on layout type (single or paired)
                    if layout_check.lower() == 'paired':
                        counts_table = counts_table.iloc[4:, [0, 3]]
                    else:
                        counts_table = counts_table.iloc[4:, [0, 2]]
                    
                    # Rename columns for clarity and append to STAR dataframes
                    counts_table.columns = ['gene', SRR]
                    star_dataframes = pd.concat([star_dataframes, counts_table], axis=1)
        
        # Remove duplicated columns if any
        star_dataframes = star_dataframes.loc[:, ~star_dataframes.columns.duplicated()]
        
        # Save combined STAR counts for the species
        star_combined_output = os.path.join(output_dir, f"combined_star_counts_{species}.csv")
        star_dataframes.to_csv(star_combined_output, index=False)
        print(f"All STAR counts for {species} combined and saved to {star_combined_output}")
        
        # Process Kallisto output
        kallisto_output_dir = os.path.join(species_output_dir, "kallisto")
        for subdir, _, files in os.walk(kallisto_output_dir):
            for file in files:
                if file == "abundance.tsv":
                    file_path = os.path.join(subdir, file)
                    
                    # Read Kallisto abundance file, which usually has columns: target_id, length, eff_length, est_counts, tpm
                    counts_table = pd.read_csv(file_path, sep='\t')
                    
                    # Extract SRA ID from the folder name
                    SRR = os.path.basename(subdir)
                    
                    # Only keep gene and estimated counts columns
                    counts_table = counts_table[['target_id', 'est_counts', 'tpm']]
                    
                    # Rename columns for clarity and append to Kallisto dataframes
                    counts = f"{SRR}_counts"
                    tpm = f"{SRR}_tpm"
                    counts_table.columns = ['gene', counts, tpm]
                    kallisto_dataframes = pd.concat([kallisto_dataframes, counts_table], axis=1)
        
        # Remove duplicated columns if any
        kallisto_dataframes = kallisto_dataframes.loc[:, ~kallisto_dataframes.columns.duplicated()]
        
        # Save combined Kallisto counts for the species
        kallisto_combined_output = os.path.join(output_dir, f"combined_kallisto_counts_{species}.csv")
        kallisto_dataframes.to_csv(kallisto_combined_output, index=False)
        print(f"All Kallisto counts for {species} combined and saved to {kallisto_combined_output}")
