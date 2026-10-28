import pandas as pd
from Bio import SeqIO
import os

def load_candidates(csv_path):
    """Load your final filtered AMP candidates."""
    df = pd.read_csv(csv_path)
    return df

def cluster_overlapping_fragments(df):
    """
    Groups fragments where one sequence is a substring of another,
    OR they share a long overlapping region (>= 15 aa common core).
    This collapses e.g. 'RSLHSIPTAFAQHSCCMPPHPCGSR' (25aa),
    'RSLHSIPTAFAQHSCCMPPHPCGSRG' (26aa), and
    'RSLHSIPTAFAQHSCCMPPHPCGSRGSR' (28aa) into ONE cluster,
    since they're all windows over the same biological region.
    """
    sequences = df['Sequence'].tolist()
    n = len(sequences)
    parent = list(range(n))  # union-find structure

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        px, py = find(x), find(y)
        if px != py:
            parent[px] = py

    # Compare every pair; union if one is a substring of the other
    for i in range(n):
        for j in range(i + 1, n):
            seq_i, seq_j = sequences[i], sequences[j]
            if seq_i in seq_j or seq_j in seq_i:
                union(i, j)

    # Assign cluster IDs
    df = df.copy()
    df['cluster_id'] = [find(i) for i in range(n)]
    return df

def pick_representatives(df):
    """
    For each cluster, keep the fragment with the HIGHEST AMP_probability
    as the representative (ties broken by longest sequence, since longer
    = more context for downstream synthesis/BLAST).
    """
    reps = (
        df.sort_values(by=['AMP_probability', 'Length'], ascending=[False, False])
          .groupby('cluster_id')
          .first()
          .reset_index()
    )
    # Also record how many overlapping fragments each cluster contained
    counts = df.groupby('cluster_id').size().rename('n_overlapping_fragments')
    reps = reps.merge(counts, on='cluster_id')
    return reps

def trace_source_protein(representatives_df, original_fasta):
    """
    For each representative candidate, find which original protein(s)
    in your immune proteome fasta it came from, by substring search.
    """
    source_map = []
    records = list(SeqIO.parse(original_fasta, "fasta"))

    for _, row in representatives_df.iterrows():
        seq = row['Sequence']
        # --- CHANGED rec.id to rec.description HERE ---
        matches = [rec.description for rec in records if seq in str(rec.seq).upper()]
        source_map.append("; ".join(matches) if matches else "NOT_FOUND")

    representatives_df = representatives_df.copy()
    representatives_df['Source_Protein'] = source_map
    return representatives_df

def main():
    input_csv = "final_results/combined_safe_amps_all_families.csv"
    original_fasta = "Labeo_rohita_immune.fasta"
    output_csv = "final_results/distinct_candidates_final.csv"

    if not os.path.exists(input_csv):
        print(f"Error: {input_csv} not found.")
        return

    df = load_candidates(input_csv)
    print(f"Loaded {len(df)} raw high-confidence fragments.")

    clustered = cluster_overlapping_fragments(df)
    n_clusters = clustered['cluster_id'].nunique()
    print(f"Collapsed into {n_clusters} distinct, non-overlapping candidate region(s).")

    representatives = pick_representatives(clustered)

    if os.path.exists(original_fasta):
        representatives = trace_source_protein(representatives, original_fasta)
    else:
        print(f"Note: {original_fasta} not found — skipping source protein trace.")

    os.makedirs("final_results", exist_ok=True)
    representatives.to_csv(output_csv, index=False)

    print(f"\n-> Saved {len(representatives)} distinct candidates to: {output_csv}\n")
    display_cols = ['Sequence', 'Length', 'AMP_family', 'AMP_probability', 
                     'n_overlapping_fragments']
    if 'Source_Protein' in representatives.columns:
        display_cols.append('Source_Protein')
    print(representatives[display_cols].to_string(index=False))

if __name__ == "__main__":
    main()