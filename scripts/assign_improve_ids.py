import pandas as pd
import argparse


def determine_delimiter(filename):
    """
    Determine the delimiter based on the file extension.
    """
    if filename.endswith('.csv'):
        return ','
    elif filename.endswith('.tsv'):
        return '\t'
    else:
        raise ValueError(f"Unsupported file type for file {filename}. Only .csv and .tsv are supported.")


def add_improve_id(previous, new, sample_col, new_name):
    """
    Add Improve IDs to a new samples.tsv file. 
    
    previous: previous sample.tsv filepath
    new: new sample.tsv filepath
    sample_col: name of column to map improve_sample_id to in the new samples.tsv
    new_name: name/save the new samples.tsv instead of overwriting it (unless you want to).
    
    Example Use:
    add_improve_id("path_to_previous_samples.tsv", "path_to_new_samples.tsv")
    """
    # Determine delimiter for the input files
    previous_delimiter = determine_delimiter(previous)
    new_delimiter = determine_delimiter(new)
    
    # Read the previous file
    previous_df = pd.read_csv(previous, sep=previous_delimiter)
    # Extract the maximum value of improve_sample_id from the previous file
    max_id = previous_df['improve_sample_id'].max()
    
    # Read the new file
    new_df = pd.read_csv(new, sep=new_delimiter)
    
    # Ceate improve_sample_id column filled with NaN values
    if 'improve_sample_id' not in new_df.columns:
        new_df['improve_sample_id'] = float('nan')
    
    # Extract unique sample_ids from the new dataframe where improve_sample_id is NaN
    unique_sample_ids = new_df[new_df['improve_sample_id'].isna()][sample_col].unique()
    
    # Create a mapping of sample_id to improve_sample_id
    id_map = {}
    for sample_id in unique_sample_ids:
        max_id += 1
        id_map[sample_id] = max_id
    
    ## Apply the mapping to the new dataframe
    new_df.loc[new_df['improve_sample_id'].isna(), 'improve_sample_id'] = new_df[sample_col].map(id_map)

    # Reorder the columns to place 'improve_sample_id' to the right of 'sample_col'
    col_idx = new_df.columns.get_loc(sample_col)
    improve_col = new_df.pop('improve_sample_id')  # Remove and get the column
    new_df.insert(col_idx + 1, 'improve_sample_id', improve_col)  # Insert it back right after sample_col
    
    # Save the updated version of the new dataframe with the same delimiter as the input 'new' file
    new_df.to_csv(new_name, sep=new_delimiter, index=False)
    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Assign improve_sample_id values to a new samples file based on a previous samples file.")
    
    parser.add_argument("-p", "--previous", type=str, required=True, help="Path to the previous samples.tsv file with existing improve_sample_id values.")
    parser.add_argument("-n", "--new", type=str, required=True, help="Path to the new samples.tsv file that needs improve_sample_id values assigned.")
    parser.add_argument("-s", "--sample_col", type=str, required=True, help="Name of column to map improve_sample_id to in the new samples.tsv.")
    parser.add_argument("-o", "--new_name", type=str, required=True, help="Name/save the new samples.tsv instead of overwriting it. You could name it the same thing if you want it overwritten though.")

    args = parser.parse_args()

    add_improve_id(args.previous, args.new, args.sample_col, args.new_name)
