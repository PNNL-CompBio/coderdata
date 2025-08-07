import gc
import polars as pl 



def main():

    # Remove Problematic Drugs before Splitting Data
    
    # Load the datasets
    all_drugs = pl.read_csv("broad_sanger_drugs.tsv", separator="\t")
    all_experiments = pl.read_csv("broad_sanger_experiments.tsv", separator="\t")

    # Define the brd_list with lowercase entries for case-insensitive matching
    brd_list = [
    'brd-k03911514',
    'brd-k07442505',
    'brd-k13185470',
    'brd-k16130065',
    'brd-k20514654',
    'brd-k27188169',
    'brd-k55473186',
    'yl54',
    'brd-k58730230',
    'brd-k79669418',
    'brd-k99584050']

    # Identify rows in all_drugs that match brd_list entries (case insensitive)
    removed_drugs = all_drugs.filter(pl.col("chem_name").str.to_lowercase().is_in(brd_list))

    # Store the improve_drug_id IDs of removed entries
    improve_drug_id = removed_drugs["improve_drug_id"].to_list()

    # Remove these rows from all_drugs and all_experiments
    all_drugs = all_drugs.filter(~pl.col("improve_drug_id").is_in(improve_drug_id))
    all_experiments = all_experiments.filter(~pl.col("improve_drug_id").is_in(improve_drug_id))
            
    all_drugs.write_csv("broad_sanger_drugs.tsv", separator="\t")
    all_experiments.write_csv("broad_sanger_experiments.tsv", separator="\t")
    
            
if __name__ == "__main__":
    main()
