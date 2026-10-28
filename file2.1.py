import pandas as pd
import glob
import os

def load_and_clean_macrel(file_path):
    """
    Loads and cleans a single Macrel prediction output file.
    """
    col_names = [
        'Access', 
        'Sequence', 
        'AMP_family', 
        'AMP_probability', 
        'Hemolytic', 
        'Hemolytic_probability'
    ]

    print(f"Loading {file_path}...")
    df = pd.read_csv(
        file_path, 
        sep='\t', 
        compression='gzip', 
        comment='#', 
        names=col_names, 
        header=None
    )
    
    # Force probability columns to be numeric
    df['AMP_probability'] = pd.to_numeric(df['AMP_probability'], errors='coerce')
    df['Hemolytic_probability'] = pd.to_numeric(df['Hemolytic_probability'], errors='coerce')
    
    # Drop any rows that became NaN (e.g., hidden headers)
    df = df.dropna(subset=['AMP_probability'])
    
    # Add sequence length column
    df['Length'] = df['Sequence'].str.len()
    
    # Reorder columns slightly so 'Length' is right next to 'Sequence'
    cols = ['Access', 'Sequence', 'Length', 'AMP_family', 'AMP_probability', 'Hemolytic', 'Hemolytic_probability']
    df = df[cols]
    
    return df

def main():
    # --- CONFIGURATION VARIABLES ---
    top_n = 15  # Change this number to show more or fewer top candidates in the terminal
    
    # 1. Find all macrel prediction files in the macrel_results_* folders
    file_pattern = "macrel_results_*/macrel.out.prediction.gz"
    file_list = glob.glob(file_pattern)
    
    if not file_list:
        print("Error: No files found matching 'macrel_results_*/macrel.out.prediction.gz'.")
        return

    # 2. Load and combine all files
    dataframes = []
    for file in sorted(file_list):
        df = load_and_clean_macrel(file)
        dataframes.append(df)
        
    # Concatenate all individual dataframes into one big dataset
    combined_df = pd.concat(dataframes, ignore_index=True)
    print(f"\nTotal combined sequences from all files: {len(combined_df)}")

    # 3. Apply Filters to the combined dataset
    highly_probable = combined_df[combined_df['AMP_probability'] >= 0.90]
    safe_amps = highly_probable[highly_probable['Hemolytic'] == 'NonHemo']
    clp_amps = safe_amps[safe_amps['AMP_family'] == 'CLP']

    print("\n--- Combined Filtering Results ---")
    print(f"Highly probable & Non-hemolytic AMPs: {len(safe_amps)}")
    print(f"Of those, Cationic Linear Peptides (CLP): {len(clp_amps)}")

    # 4. Create the final results directory
    output_dir = "final_results"
    os.makedirs(output_dir, exist_ok=True)
    print(f"\nCreated folder: '{output_dir}/'")

    # 5. Save the combined results to CSV
    safe_amps_path = os.path.join(output_dir, "combined_safe_amps.csv")
    safe_amps.to_csv(safe_amps_path, index=False)
    print(f"-> Saved: {safe_amps_path}")

    clp_amps_path = os.path.join(output_dir, "combined_clp_amps.csv")
    clp_amps.to_csv(clp_amps_path, index=False)
    print(f"-> Saved: {clp_amps_path}")

    # 6. Display top N results across all batches
    if not safe_amps.empty:
        print(f"\nTop {top_n} non-hemolytic AMP candidates overall (sorted by AMP probability):")
        # Uses the top_n variable dynamically here
        top_candidates = safe_amps.sort_values(by='AMP_probability', ascending=False).head(top_n)
        print(top_candidates[['Access', 'Sequence', 'Length', 'AMP_family', 'AMP_probability']])

if __name__ == "__main__":
    main()