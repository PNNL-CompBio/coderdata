import pandas as pd

# Read proteomics data
p_data = pd.read_csv("proteomics.csv.gz")

# Read samples data
s_data = pd.read_csv("samples.csv")

# Read genes data
g_data = pd.read_csv("genes.csv")
g_data.rename(columns={"other_id": "ENSG_col"}, inplace=True)
g_data = g_data.groupby('ENSG_col').first().reset_index()
g_data = g_data.dropna(subset=['ENSG_col'])


# Merging operations
merged = p_data.merge(s_data[['improve_sample_id', 'other_id']], on='improve_sample_id', how='left')
merged = merged.merge(g_data[["entrez_id", "gene_symbol", "ENSG_col"]], on="entrez_id", how="left")

agg_df = merged.groupby(['improve_sample_id', 'entrez_id', "gene_symbol", "ENSG_col"]).agg({
    'proteomics': 'mean',
    'other_id': 'first'
}).reset_index()

agg_df.drop(columns=["improve_sample_id"], inplace=True)
agg_df.drop_duplicates(inplace=True)
agg_df = agg_df.dropna(subset=['other_id'])

# Pivot merged dataframe
pivoted = agg_df.pivot(index='other_id', columns=["ENSG_col", 'entrez_id', "gene_symbol"], values='proteomics').reset_index()

cols = list(pivoted.columns)
cols[0] = ('', '', '')
header = pd.MultiIndex.from_tuples(cols)
pivoted.columns = header

# Transpose and reset multi-index
pivoted = pivoted.T
pivoted = pivoted.reset_index().drop_duplicates(subset=['level_1', 'level_2']).set_index(['level_0', 'level_1', 'level_2'])
pivoted = pivoted.T

# Final adjustment to the column headers
cols = list(pivoted.columns)
cols[0] = ('', '', '')
header = pd.MultiIndex.from_tuples(cols)
pivoted.columns = header

# Save to CSV
pivoted.to_csv("proteomics_restructure.tsv", index=False, sep="\t")
