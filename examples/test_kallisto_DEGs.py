import pandas as pd
import os
#%%
species_output_dir = "/home/mpimp-golm.mpg.de/khayer6872/AG-Gutjahr/Group_Projects/abul_k/P072024_ZH_AK_RNAseq/outputs/kartikey_rnaseq"
kallisto_output_dir = species_output_dir
kallisto_dataframes = pd.DataFrame()
for subdir, _, files in os.walk(kallisto_output_dir):
    for file in files:
        if file.endswith(".tsv"):
            file_path = os.path.join(subdir, file)
            
            # Read Kallisto abundance file, which usually has columns: target_id, length, eff_length, est_counts, tpm
            counts_table = pd.read_csv(file_path, sep='\t')
            
            # Extract SRA ID from the folder name
            SRR = os.path.basename(file_path)
            SRR = SRR.split('.')[0]
            # Only keep gene and estimated counts columns
            counts_table = counts_table[['target_id', 'est_counts', 'tpm']]
            
            # Rename columns for clarity and append to Kallisto dataframes
            counts = f"{SRR}_counts"
            tpm = f"{SRR}_tpm"
            counts_table.columns = ['gene', counts, tpm]
            kallisto_dataframes = pd.concat([kallisto_dataframes, counts_table], axis=1)
#%%
# Remove duplicated columns if any
kallisto_dataframes = kallisto_dataframes.loc[:, ~kallisto_dataframes.columns.duplicated()]
#%%
