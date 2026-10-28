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
    
    # Reorder columns slightly
    cols = ['Access', 'Sequence', 'Length', 'AMP_family', 'AMP_probability', 'Hemolytic', 'Hemolytic_probability']
    df = df[cols]
    
    return df

def main():
    # --- CONFIGURATION VARIABLES ---
    top_n = 15             # Number of top candidates to display in the terminal
    min_probability = 0.90 # Minimum AMP probability threshold
    
    # 1. Find all macrel prediction files
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
        
    combined_df = pd.concat(dataframes, ignore_index=True)
    print(f"\nTotal combined sequences from all files: {len(combined_df)}")

    # 3. Apply Filters
    highly_probable = combined_df[combined_df['AMP_probability'] >= min_probability]
    safe_amps = highly_probable[highly_probable['Hemolytic'] == 'NonHemo']

    print("\n--- Combined Filtering Results ---")
    print(f"Highly probable (>= {min_probability}) & Non-hemolytic AMPs: {len(safe_amps)}")

    # 4. Create the final results directory
    output_dir = "final_results"
    os.makedirs(output_dir, exist_ok=True)
    print(f"\nCreated folder: '{output_dir}/'")

    # 5. Save the master combined results
    safe_amps_path = os.path.join(output_dir, "combined_safe_amps_all_families.csv")
    safe_amps.to_csv(safe_amps_path, index=False)
    print(f"-> Saved Master File: {safe_amps_path}")

    # 6. DYNAMICALLY save separate CSVs for each AMP family found
    if not safe_amps.empty:
        # Get a list of all unique families present in the filtered data
        found_families = safe_amps['AMP_family'].unique()
        print(f"\nFound {len(found_families)} distinct AMP family type(s): {', '.join(found_families)}")
        
        for family in found_families:
            # Filter the dataframe for just this family
            family_df = safe_amps[safe_amps['AMP_family'] == family]
            
            # Create a clean filename (e.g., combined_cdp_amps.csv)
            family_filename = f"combined_{family.lower()}_amps.csv"
            family_path = os.path.join(output_dir, family_filename)
            
            # Save it
            family_df.to_csv(family_path, index=False)
            print(f"-> Saved {family} family ({len(family_df)} sequences): {family_path}")

        # 7. Display top N results
        print(f"\nTop {top_n} non-hemolytic AMP candidates overall:")
        top_candidates = safe_amps.sort_values(by='AMP_probability', ascending=False).head(top_n)
        print(top_candidates[['Access', 'Sequence', 'Length', 'AMP_family', 'AMP_probability']])
    else:
        print("\nNo sequences passed the current probability and non-hemolytic filters.")

if __name__ == "__main__":
    main()


    