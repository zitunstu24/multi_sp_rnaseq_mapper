#%%
import os
import pandas as pd
#%%
df = pd.read_csv("/winmounts/khayer6872/AG-Gutjahr/Group_Projects/abul_k/playground/sunflower_AM_AF_v1/inputs/master_table.csv")
output_dir = "/winmounts/khayer6872/AG-Gutjahr/Group_Projects/abul_k/playground/sunflower_AM_AF_v1/outputs/tmp"

#%%
for i in df['species'].unique():
        
        print(f"started combining for {i}")

        species_output_dir = os.path.join(output_dir, i)

        all_dataframes = pd.DataFrame()

        for subdir, _, files in os.walk(species_output_dir):

            for file in files:

                if file.endswith("ReadsPerGene.out.tab"):
                    file_path = os.path.join(subdir, file)

                    counts_table = pd.read_csv(file_path,header=0,delimiter="\t")

                    SRR = os.path.basename(subdir)

                    print(SRR)

                    layout_check = df[df['SRA_ID'] == SRR]['layout']

                    layout_check = layout_check.iloc[0] if isinstance(layout_check, pd.Series) else layout_check

                    if layout_check.lower() == 'paired':
                        counts_table = counts_table.iloc[4:, [0,3]]
                        counts_table = counts_table.rename(
                            columns={
                                    counts_table.columns[0] : 'gene',
                                    counts_table.columns[1] : SRR
                                }
                        )

                    else:
                        counts_table = counts_table.iloc[4:, [0,2]]
                        counts_table = counts_table.rename(
                            columns={
                                    counts_table.columns[0] : 'gene',
                                    counts_table.columns[1] : SRR
                                }
                        )

                    print(counts_table.head(2))

                    print(type(counts_table))

                    all_dataframes = pd.concat([all_dataframes,counts_table], axis=1)
                    

                    # print(all_dataframes.head(2))
                    #all_dataframes.merge(counts_table, how='cross')
        print(f'finished for {i}')
        
        print(all_dataframes.head(2))

        all_dataframes = all_dataframes.loc[:, ~all_dataframes.columns.duplicated()]

        all_dataframes.to_csv(os.path.join(output_dir, f"combined_counts_{i}.csv"), index=False)

# %%
