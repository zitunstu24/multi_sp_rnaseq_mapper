#%%
from multi_sp_rnaseq_mapper.config import Config
import os
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
#%%
output_dir = '/home/mpimp-golm.mpg.de/khayer6872/p0424_AK_KAR/outputs'
# master_table = config.get_master_table()
df = pd.read_csv('/home/mpimp-golm.mpg.de/khayer6872/p0424_AK_KAR/inputs/master_table.csv')
df = df[df['species'] == 'marchantia']


#%%

for i in df['species'].unique():

    print(f"started DEGs for {i}")

    filtered_df = df[df['species'] == i]

    counts_df_path = os.path.join(output_dir, f"combined_counts_{i}.csv")
    counts_df = pd.read_csv(counts_df_path)

    all_dataframes = pd.DataFrame()

    for sample_groups in filtered_df['groups'].unique():
        
        group_controls = filtered_df[filtered_df['groups'] == sample_groups]['respective_controls']

        if group_controls.isna().all() or (group_controls == '').all():

            print(f"Group '{sample_groups}' is controls group. Moving to the next treatment group...")

            continue

        else:
            print(f"Group '{sample_groups}' has valid controls. comparing it with {group_controls}...")

            
            group_controls = group_controls.tolist()
            filtered_subset = filtered_df[(filtered_df['groups'] == sample_groups) |
                            (filtered_df['groups'] == group_controls[0])]
            
            filtered_subset['custom_sort'] = filtered_subset['groups'].apply(lambda x: 0 if x == sample_groups else 1)
            filtered_subset = filtered_subset.sort_values(by='custom_sort')
            filtered_subset.drop(columns=['custom_sort'], inplace=True)
            filtered_subset.set_index('SRA_ID', inplace=True)
            filtered_subset.index.name = None

        SRR_treatment = filtered_df[filtered_df['groups'] == sample_groups]['SRA_ID'].tolist()
        SRR_control = filtered_df[filtered_df['groups'] == group_controls[0]]['SRA_ID'].tolist()
        columns_to_keep = ['gene'] + SRR_treatment + SRR_control



        counts_subset = counts_df.loc[:, columns_to_keep]

        counts_subset_processed = counts_subset.set_index('gene')
        counts_subset_processed.index.name = None
        counts_subset_processed = counts_subset_processed.transpose() 
        counts_subset_processed.fillna(0, inplace=True)
        counts_subset_processed = counts_subset_processed.reindex(filtered_subset.index)

        inference = DefaultInference(n_cpus=100)
        dds = DeseqDataSet(
                counts=counts_subset_processed,
                metadata=filtered_subset,
                design_factors="groups",
                refit_cooks=True,
                inference=inference,
            )

        dds.deseq2()


        stat_res = DeseqStats(dds, inference=inference)
        stat_res.summary()

        result = stat_res.results_df
        result['gene'] = result.index

        suffix = "_" + str(group_controls[0]) + "_vs_" + str(sample_groups)

        result = result.add_suffix(suffix)

        all_dataframes = pd.concat([all_dataframes,result],axis = 1)

    print(f"done DEGs for {i}")
    gene_columns = [col for col in all_dataframes.columns if 'gene' in col]

    gene_column_to_keep = gene_columns[0]

    all_dataframes = all_dataframes[[gene_column_to_keep] + [col for col in all_dataframes.columns if col not in gene_columns]]

    cols = [gene_column_to_keep] + [col for col in all_dataframes.columns if col != gene_column_to_keep]
    all_dataframes = all_dataframes[cols]
    all_dataframes = all_dataframes.rename(columns={gene_column_to_keep: 'gene'})

    all_dataframes.to_csv(os.path.join(output_dir, f"test_combined_DEGs_{i}.csv"), index=False)
