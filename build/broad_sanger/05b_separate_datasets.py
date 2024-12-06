import gc
import polars as pl 



def main():
    datasets_to_process = ["CCLE", "CTRPv2", "PRISM", "GDSCv1", "GDSCv2", "FIMM", "gCSI", "NCI60"]
    omics_datatypes = ["transcriptomics","proteomics", "copy_number","mutations"] # csv 
    samples_datatypes = ["samples"] #csv

    drugs_datatypes = ["drugs", "drug_descriptors"] # tsv


    dataset_sources = {
        "CCLE": ["Broad"],
        "CTRPv2": ["Broad"],
        "PRISM": ["Broad"],
        "GDSCv1": ["Sanger"],
        "GDSCv2": ["Sanger"],
        "FIMM": ["Broad"],
        "gCSI": ["Broad"],  # gCSI generates its own omics data but it is comparable to CCLE. In future, retrive gCSI omics.
        "NCI60": ["Broad"]
    }

    for dataset in datasets_to_process:
        exp = pl.read_csv("broad_sanger_experiments.tsv", separator="\t") # Keeping memory down, so I will not be making copies.
        exp = exp.filter(pl.col("study") == dataset)

        # Extract information to separate out datasets
        exp_improve_sample_ids = exp["improve_sample_id"].unique().to_list()
        exp_improve_drug_ids = exp["improve_drug_id"].unique().to_list()

        # Write Filtered Experiments File to TSV. Then delete it from memory.
        exp_filename = f"/tmp/{dataset}_experiments.tsv".lower()
        exp.write_csv(exp_filename, separator="\t")
        del exp
        gc.collect()


        #Filter Samples files, write to file, delete from mem.
        for samples in samples_datatypes:
            samples_filename_in = f"broad_sanger_{samples}.csv"
            samples_filename_out = f"/tmp/{dataset}_{samples}.csv".lower()
            samples_df = pl.read_csv(samples_filename_in)
            samples_df = samples_df.filter(pl.col("improve_sample_id").is_in(exp_improve_sample_ids))
            samples_df.write_csv(samples_filename_out) #csv
            del samples_df
            gc.collect()

        #One by one, filter other Omics files, write to file, delete from mem.
        for omics in omics_datatypes:
            omics_filename_in = f"broad_sanger_{omics}.csv"
            omics_filename_out = f"/tmp/{dataset}_{omics}.csv".lower()
            omics_df = pl.read_csv(omics_filename_in)
            omics_df = omics_df.filter(pl.col("improve_sample_id").is_in(exp_improve_sample_ids))
            omics_df = omics_df.filter(pl.col("source").is_in(dataset_sources[dataset]))
            omics_df.write_csv(omics_filename_out) #csv
            del omics_df
            gc.collect()


        #One by one, filter other Drugs files, write to file, delete from mem.
        for drugs in drugs_datatypes:
            drugs_filename_in = f"broad_sanger_{drugs}.tsv"
            drugs_filename_out = f"/tmp/{dataset}_{drugs}.tsv".lower()
            if drugs == "drug_descriptors":
                drugs_df = pl.read_csv(drugs_filename_in,separator="\t",
                                       dtypes={"improve_drug_id": pl.Utf8,
                                                         "structural_descriptor": pl.Utf8,
                                                         "descriptor_value": pl.Utf8}
                                      )

            else:
                drugs_df = pl.read_csv(drugs_filename_in,separator="\t")

            drugs_df = drugs_df.filter(pl.col("improve_drug_id").is_in(exp_improve_drug_ids))
            drugs_df.write_csv(drugs_filename_out,separator="\t") #tsv
            del drugs_df
            gc.collect()
            
if __name__ == "__main__":
    main()
