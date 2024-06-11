# Imports:
import pandas as pd
import argparse

# Known modification positions and gene differences for 16s and 23s rRNA. 
known_mod_positions23s = {
    747: "m¹G", 748: "Ψ", 749: "m⁵U", 957: "Ψ", 1620: "m⁶A",
    1837: "m²G", 1915: "Ψ", 1919: "m³Ψ", 1921: "Ψ", 1943: "m⁵U",
    1966: "m⁵C", 2034: "m⁶A", 2073: "m⁷G", 2255: "Gm", 2449: "m²G",
    2453: "D", 2461: "Ψ", 2502: "Cm", 2505: "OH⁵C", 2507: "m²A",
    2508: "Ψ", 2556: "Um", 2584: "Ψ", 2608: "Ψ", 2609: "Ψ"
}

known_mod_positions16s = {
    516: "Ψ", 527: "m⁷G", 966: "m²G", 967: "m⁵C", 1207: "m²G",
    1402: "m⁴Cm", 1407: "m⁵C", 1498: "m³U", 1518: "m⁶₂A", 1519: "m⁶₂A"
}

# Define the gene difference positions for 23S and 16S
genediff_positions23s = [83, 137, 142, 264, 356, 357, 388, 541, 542, 543, 544, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 1213, 1476, 1725, 1727, 1728, 1729, 1732, 1735, 1736, 1737, 1867, 1872, 1873, 1874, 1876, 2136, 2167, 2207, 2209, 2213, 2219, 2798, 2800, 2803, 2806, 2907, 2908, 2909]
genediff_positions16s = [79, 90, 210, 449, 1002, 1006, 1010, 1019, 1020, 1021, 1022, 1023, 1038]



def process_and_combine_EpinNano_csv(file1_path, file2_path, output_file):
    """
    Process and combine two EpinNano CSV files into a single CSV file.

    Args:
        file1_path (str): Path to the first input CSV file.
        file2_path (str): Path to the second input CSV file.
        output_file (str): Path to the output CSV file.

    Returns:
        None, but saves the combined and sorted CSV file to the specified path and with the specified name.

    Raises:
        None

    This function loads two EpinNano Differ.R Output CSV files, determines the structure of the files,
    renames columns for uniformity, merges the dataframes, expands the 'chr_pos' column,
    fills missing columns with default values, sorts the dataframe by 'position', and
    saves the final dataframe to a CSV file.

    Note: The 'chr_pos' format is assumed to be like "ReferenceRNA 100 T +". Adjust the regular
    expression if your data varies.

    """
    # Load the files
    df1 = pd.read_csv(file1_path)
    df2 = pd.read_csv(file2_path)
    
    # Determine which file is which
    if 'delta_sum_err' in df1.columns:
        delta_df = df1
        linear_regression_df = df2
    else:
        delta_df = df2
        linear_regression_df = df1
    
    # Rename columns for uniformity
    rename_columns = {
        'ko_feature': 'ko_sum_err', 
        'wt_feature': 'wt_sum_err'
    }
    delta_df.rename(columns=rename_columns, inplace=True)
    linear_regression_df.rename(columns=rename_columns, inplace=True)
    
    # Merge dataframes
    merged_df = pd.merge(delta_df, linear_regression_df, on=['chr_pos', 'ko_sum_err', 'wt_sum_err'])
    
    # Correct regular expression to match 'chr_pos' format accurately
    # Assuming format is like "23rRNA 100 T +", adjust if your data varies
    chr_pos_expanded = merged_df['chr_pos'].str.extract(r'(\d+([sr])RNA) (\d+) ([ATCGN]) ([+-])')
    
    if chr_pos_expanded.isnull().values.any():
        print("Warning: Some rows could not be split correctly. Please check the 'chr_pos' format.")
    
    chr_pos_expanded.columns = ['rRNA', 'position', 'base', 'direction']
    
    # Convert 'position' to numeric
    chr_pos_expanded['position'] = pd.to_numeric(chr_pos_expanded['position'])
    
    # Concatenate expanded 'chr_pos' with the rest of the DataFrame
    result_df = pd.concat([chr_pos_expanded, merged_df.drop(columns=['chr_pos'])], axis=1)
    
    # Define the final columns order
    final_columns = ['rRNA', 'position', 'base', 'direction', 'ko_sum_err', 'wt_sum_err', 'delta_sum_err', 'z_scores', 'z_score_prediction', 'lm_Bonferroni_outlier_test', 'lm_residuals', 'lm_residuals_z_score', 'lm_residuals_z_scores_prediction']
    
    # Ensure all required columns are present, fill missing with default value if necessary
    for col in final_columns:
        if col not in result_df.columns:
            result_df[col] = None  # or a suitable default value
    
    # Sort by 'position'
    final_df = result_df.sort_values(by='position')
    
    # Save to CSV
    final_df[final_columns].to_csv(output_file, index=False)
    
    print(f"Combined and sorted CSV saved as {output_file}")

# Example usage
# process_and_combine_EpinNano_csv('path/to/delta-sum_err.prediction.csv', 'path/to/linear-regression.prediction.csv', 'path/to/combined_output.csv')



def label_positions(row, known_mods, gene_diffs):
    """
    Labels the positions based on their characteristics.
    Use 5-mer window around the position to determine the label.
    
    Args:
        row (pandas.Series): The row containing the position to be labeled.
        known_mods (list): List of known modifications.
        gene_diffs (list): List of gene differences.

    Returns:
        str: The label for the position.
    """
    pos = row['position']
    if pos in known_mods:
        return 'Known Modification'
    elif pos in gene_diffs:
        return 'Gene Difference'
    elif any(pos - 2 <= diff_pos <= pos + 2 for diff_pos in gene_diffs) and any(pos - 2 <= mod_pos <= pos + 2 for mod_pos in known_mods):
        return 'Near Mod and Diff'
    elif any(pos - 2 <= diff_pos <= pos + 2 for diff_pos in gene_diffs):
        return 'Near Diff'
    elif any(pos - 2 <= mod_pos <= pos + 2 for mod_pos in known_mods):
        return 'Near Mod'
    else:
        return 'Canonical Base'
    


def add_mod_id(row, known_mods):
    """
    Adds a modification ID to the given row based on the known modifications dictionary.

    Args:
        row (dict): The row containing the position of the modification.
        known_mods (dict): A dictionary mapping positions to modification IDs.

    Returns:
        str: The modification ID corresponding to the position in the row, or an empty string if not found.
    """
    return known_mods.get(row['position'], '')



def process_Epinano_rRNA(input_file):
    """
    Process the Epinano rRNA data from the given processed CSV input file.
    
    Args:
        input_file (str): The path to the input CSV file containing the Epinano rRNA data.
        
    Returns:
        A new CSV file containing the processed data for modifications and gene differences.
        
    Raises:
        None
    """
    df = pd.read_csv(input_file)
    
    # Determine and apply labels
    df['pos_id'] = df.apply(lambda row: label_positions(row, known_mod_positions16s if row['rRNA'] == '16sRNA' 
                                                        else known_mod_positions23s, genediff_positions16s 
                                                        if row['rRNA'] == '16sRNA' 
                                                        else genediff_positions23s), axis=1)
    df['mod_id'] = df.apply(lambda row: add_mod_id(row, known_mod_positions16s if row['rRNA'] == '16sRNA' 
                                                   else known_mod_positions23s), axis=1)
    
    # Save to new CSV
    output_file = input_file.replace('.csv', '_Processed.csv')
    df.to_csv(output_file, index=False)
    print(f"Processed file saved as {output_file}")

# Example usage
# process_Epinano_rRNA('combined_output.csv')
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process Epinano replicate files and combine or label them for further analysis.")
    parser.add_argument('--mode', type=str, choices=['combine', 'label'], required=True, help="Mode of operation: 'combine' to merge and process files, 'label' to add labels based on known modifications and differences.")
    parser.add_argument('--file1', type=str, help="Path to the first input CSV file for combining. Required if mode is 'combine'.")
    parser.add_argument('--file2', type=str, help="Path to the second input CSV file for combining. Required if mode is 'combine'.")
    parser.add_argument('--output', type=str, help="Path to the output CSV file. Required if mode is 'combine'.")
    parser.add_argument('--input', type=str, help="Path to the input CSV file for labeling. Required if mode is 'label'.")
    args = parser