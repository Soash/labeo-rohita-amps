# In Silico Mining of Antimicrobial Peptide Candidates from the Immune Proteome of *Labeo rohita* Reveals a Cryptic Transferrin-Derived Peptide

**Author(s):** Md Abu Sadat Soash  
*Department of Genetic Engineering and Biotechnology*<br>
*Shahjalal University of Science and Technology, Sylhet, Bangladesh*  
**Correspondence:** [mdabu13@student.sust.edu](mailto:mdabu13@student.sust.edu)

[![UniProt Proteome](https://img.shields.io/badge/UniProt-Labeo%20rohita%20(BAU--BD--2019)-005A9C.svg)](https://www.uniprot.org/taxonomy/84645)
[![Macrel](https://img.shields.io/badge/Classifier-Macrel%20v1.4.5-teal.svg)](https://github.com/BigDataBiology/macrel)
[![SeqKit](https://img.shields.io/badge/FASTA%20Toolkit-SeqKit%20v2.0+-orange.svg)](https://bioinf.shenwei.me/seqkit/)
[![Python](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

---

## Abstract

* **Objective:** Antimicrobial peptides (AMPs) are short, amphipathic molecules occurring in innate immune systems across all domains of life that destabilize microbial membranes and modulate host immune responses. Increasingly, active AMPs are identified as "cryptic" peptides encrypted within larger functional host proteins. In South Asian aquaculture, the escalating crisis of antimicrobial resistance (AMR) severely compromises production of *Labeo rohita* (rohu), yet host-derived AMPs remain largely unexplored in this high-value aquaculture species. This study aims to systematically discover and characterize candidate AMP sequences encrypted within the immune proteome of *Labeo rohita*.
* **Methodology:** The reference proteome for *Labeo rohita* strain BAU-BD-2019 (specimen collected from Rangpur, Bangladesh; UniProt TaxID: 84645) was filtered for proteins annotated with immune system processes ($n = 606$). An exhaustive sliding-window fragmentation algorithm generated all possible $15\text{--}50$ amino acid subsequences ($10,049,499$ unique fragments). The library was partitioned using **SeqKit** and screened across 28 parallel threads using **Macrel**, a 22-descriptor random forest AMP and hemolytic activity classifier. Predictions were filtered for high confidence ($\text{AMP Probability} \ge 0.90$) and predicted non-hemolytic safety ($\text{Hemolytic} = \text{'NonHemo'}$). Redundant sliding-window overlaps were collapsed using a union-find graph clustering algorithm and traced back to their source proteins.
* **Findings:** Initial proteome-wide screening classified **251,328 fragments** as AMP-positive. Applying stringent confidence thresholds and non-hemolytic filters condensed the pool to 9 high-confidence hits across all 3 partitions, collapsing into **four distinct candidate regions**. The single most physiologically viable candidate is a **40-residue cationic, disulfide-stabilized peptide** derived from serotransferrin (`rohu_fragment_503000`), carrying a net charge of $+4.0$, high aqueous solubility, an amphiphilic index of $0.86$, and an AMP probability of $0.901$ with low predicted hemolytic liability ($0.347$). Canonical basic protease cleavage sites flank the peptide, consistent with endogenous liberation by inflammatory proteases. Sequence homology comparison against curated peptide repositories identified two validated antimicrobial homologs, **sOT2** and **cOT2** (from softshell turtle and Siamese crocodile ovotransferrin), sharing the conserved `KSCHTGL` active core motif. The three remaining candidate regions (including a GATA-2 zinc-finger domain and an E3 ligase) reflect classifier artifacts driven by high cysteine and basic residue density rather than true secretable effector peptides.
* **Implications:** Proteome-wide machine learning screening requires rigorous source-protein contextualization to filter out structural false positives. The discovered cryptic serotransferrin-derived peptide constitutes a high-priority lead for peptide synthesis, structural characterization, and downstream formulation of species-relevant antimicrobials for sustainable aquaculture.

---

## 1. Introduction & Biological Background

Antimicrobial resistance (AMR) is an acute threat to aquaculture sustainability in South Asia, where intensive fish farming relies heavily on antibiotics to suppress outbreaks of opportunistic pathogens such as *Aeromonas hydrophila*, *Edwardsiella tarda*, and *Pseudomonas fluorescens*. The widespread emergence of multidrug-resistant bacteria poses severe risks to food safety, economic security, and public health.

Antimicrobial peptides (AMPs) represent an attractive therapeutic paradigm due to their rapid, membrane-disruptive mechanisms of action that minimize the likelihood of bacterial resistance evolution. While classical AMPs are expressed as dedicated precursor prepropeptides (e.g., defensins, cathelicidins, hepcidins), growing evidence indicates that innate immune defenses also exploit **cryptic antimicrobial peptides**—bioactive sequences embedded within larger functional proteins (such as lactoferrin, ovotransferrin, and hemoglobin) that are released upon targeted proteolytic cleavage during tissue injury or inflammation.

```
          Aquaculture AMR Crisis in South Asia (Labeo rohita)
                                   │
                   ┌───────────────┴───────────────┐
                   ▼                               ▼
       Excessive Antibiotic Use &         Vulnerability to Opportunistic
      Emergence of Resistant Strains        Pathogens (Aeromonas, etc.)
                   │                               │
                   └───────────────┬───────────────┘
                                   ▼
          In Silico Mining of Cryptic Innate Host Defense Peptides
                                   │
                   ┌───────────────┴───────────────┐
                   ▼                               ▼
       606 Immune Proteins Curated        10,049,499 Overlapping Fragments
        from BAU-BD-2019 Proteome           (15–50 aa Sliding Window)
                   │                               │
                   └───────────────┬───────────────┘
                                   ▼
             Machine Learning Screening via Macrel (28 Threads)
                   │
                   ▼
         251,328 AMP-Positive Predictions Screened
                   │
                   ▼
         Stringent Filtering (Probability ≥ 0.90 & Non-Hemolytic)
                   │
                   ▼
         Union-Find Clustering (4 Distinct Non-Overlapping Regions)
                   │
                   ▼
         Discovery of Cryptic Transferrin-Derived Bioactive AMP
```

1. **Host-Derived Cryptic Peptides:** Unlike exogenous synthetic antimicrobials, cryptic peptides derived from endogenous host proteins exhibit low immunogenicity, negligible hemolytic toxicity, and species-specific compatibility.
2. **Transferrin Family Antimicrobial Crypticity:** Transferrins and lactoferrins are primary iron-chelating glycoproteins of the innate immune system. In vertebrates, proteolytic cleavage of the N-terminal transferrin domain releases potent microbicidal peptides (such as lactoferricin in mammals, and sOT2/cOT2 in reptiles).
3. **Study Aims:** This study establishes an exhaustive computational screening framework to survey the entire annotated immune proteome of *Labeo rohita* (rohu), isolate candidate cryptic peptides using random forest models, and validate their biological authenticity against structural and evolutionary benchmarks.

---

## 2. Dataset & Target Proteome Stratification

The reference proteome for *Labeo rohita* strain **BAU-BD-2019** (specimen collected from Rangpur, Bangladesh; UniProt Taxonomy ID: [84645](https://www.uniprot.org/taxonomy/84645)) was retrieved from UniProt and subset for immune system processes:

| Dataset Metric | Verified Statistic | Methodological / Analytical Significance |
| :--- | :---: | :--- |
| **Reference Organism** | *Labeo rohita* (Rohu) | Major aquaculture carplet species in South Asia |
| **Proteome Strain / Origin** | **BAU-BD-2019** (Rangpur) | High-quality reference proteome assembly |
| **Curated Immune Proteins** | **606** | UniProt entries annotated for immune system processes |
| **Immune Protein Length Range** | **51 – 5,653 aa** | Mean length: $511.4$ amino acids |
| **Sliding Window Range** | **15 – 50 aa** | Exhaustive single-residue step fragmentation |
| **Total Unique Peptide Fragments** | **10,049,499** | Exhaustive candidate fragment space |
| **Parallelization Splits** | **3 Partitions** | $\sim 3,349,833$ sequences per partition via SeqKit |
| **Macrel AMP-Positive Predictions** | **251,328** | Total raw hits flagged as potential AMPs |
| **High-Confidence Filter** | **AMP Prob $\ge$ 0.90** | Random forest ensemble classification threshold |
| **Biosafety Filter** | **Non-Hemolytic** | $\text{Hemolytic} = \text{'NonHemo'}$, probability $< 0.50$ |
| **Filtered Safe Peptides** | **9 fragments** | Sequences meeting joint probability and safety criteria |
| **Distinct Candidate Clusters** | **4 regions** | Collapsed non-overlapping biological regions |

---

## 3. In Silico Screening & Computational Pipeline

```
[ Labeo_rohita_immune.fasta ] (606 proteins)
               │
               ▼  file1.py (Sliding Window: 15–50 aa)
[ Labeo_rohita_AMPs.fasta ] (10,049,499 unique fragments, ~540 MB)
               │
               ▼  seqkit split2 -p 3
[ splits/part_001.fasta ]   [ splits/part_002.fasta ]   [ splits/part_003.fasta ]
               │                           │                           │
               ▼                           ▼                           ▼
      Macrel (28 threads)         Macrel (28 threads)         Macrel (28 threads)
               │                           │                           │
               ▼                           ▼                           ▼
[ macrel_results_1/ ]       [ macrel_results_2/ ]       [ macrel_results_3/ ]
               │                           │                           │
               └───────────────────────────┼───────────────────────────┘
                                           │
                                           ▼  file2.1.py / file2.2.py
                       [ final_results/combined_safe_amps_all_families.csv ]
                                           │
                                           ▼  dedup_candidates.py (Union-Find)
                       [ final_results/distinct_candidates_final.csv ]
                                           │
                                           ▼
                     Structural Modeling & Evolutionary Validation
```

### 3.1 Phase 1: In Silico Proteolytic Sliding Window Fragmentation (`file1.py`)
Each immune-related protein sequence was processed with an exhaustive sliding window spanning lengths $L \in [15, 50]$. All canonical 20-amino-acid fragments were deduplicated in memory and indexed as `rohu_fragment_X`, producing $10,049,499$ unique fragments.

### 3.2 Phase 2: High-Throughput Classification via Macrel
Due to the massive sequence space, the peptide library was partitioned into 3 equal batches using `seqkit split2` and screened with **Macrel v1.4.5** across 28 hardware threads:
* Macrel calculates 22 physicochemical and structural features including charge, hydrophobic moment, isoelectric point, and secondary structure propensity.
* Features are evaluated through an ensemble of 22 random forest classifiers trained on validated antimicrobial and hemolytic peptide datasets.

### 3.3 Phase 3: Toxicological and Stringency Filtering (`file2.1.py` / `file2.2.py`)
Macrel output tables (`macrel.out.prediction.gz`) were parsed and filtered dynamically:
1. **AMP Probability:** Set strictly to $\ge 0.90$.
2. **Hemolytic Safety:** Restricted to `NonHemo` (predicted hemolytic probability $< 0.50$).
3. **Family Stratification:** Filtered peptides were partitioned into mechanistic families (CLP: Cationic Linear Peptides; CDP: Cationic Disulfide-containing Peptides). All top-tier hits belonged to the CDP class.

### 3.4 Phase 4: Union-Find Graph Clustering (`dedup_candidates.py`)
Because sliding windows generate overlapping subsequences over the same biological locus (e.g., 25-mer, 26-mer, and 28-mer variants of the same epitope), a union-find disjoint-set algorithm was implemented:
* Sequences sharing a substring relation or common core region ($\ge 15$ aa) were collapsed into a single cluster.
* The representative peptide for each cluster was selected by maximizing $\text{AMP Probability}$ (ties resolved by sequence length).

### 3.5 Phase 5: Reverse Source-Protein Tracing & Contextualization
Each representative cluster sequence was mapped back to its parent protein in `Labeo_rohita_immune.fasta` to establish subcellular localization, domain architecture, and physiological cleavage feasibility.

---

## 4. Screening Results & Candidate Characterization

### 4.1 Four Distinct Candidate Regions

The pipeline resolved the $10,049,499$ fragments into **four distinct candidate regions**:

| Cluster ID | Accession ID | Representative Sequence | Length | AMP Class | AMP Prob | Hemolytic Status (Prob) | Overlapping Windows | Source Protein Context |
| :---: | :--- | :--- | :---: | :---: | :---: | :---: | :---: | :--- |
| **2** | `rohu_fragment_503000` | `KSCHTGLNRNAGWKVPDSAICGKTPDCTLYKFFSKGCAPG` | **40 aa** | **CDP** | **0.901** | NonHemo (0.347) | 1 | **Serotransferrin** (`A0ABQ8MX81_LABRO`) |
| **4** | `rohu_fragment_3727390` | `RSLHSIPTAFAQHSCCMPPHPCGSR` | **25 aa** | **CDP** | **0.960** | NonHemo (0.455) | 4 | **E3 ubiquitin ligase TRAF3IP2** (`A0ABQ8MLE2_LABRO`) |
| **6** | `rohu_fragment_5729201` | `DGYRSQKSFCKIFGRSLCPGTSCGA` | **25 aa** | **CDP** | **0.931** | NonHemo (0.376) | 1 | **T-cell surface glycoprotein CD1a** (`A0ABQ8L5R9_LABRO`) |
| **7** | `rohu_fragment_3831189` | `KSRAFSGRECVNCGATSTPLWRRDGTGHYLCNACGLYHKMNGQ` | **43 aa** | **CDP** | **0.901** | NonHemo (0.317) | 2 | **GATA-binding factor 2** (`A0ABQ8MMZ7_LABRO`) |

---

### 4.2 Priority Candidate 1: Cryptic Serotransferrin-Derived Bioactive Peptide

The 40-residue peptide from *Labeo rohita* serotransferrin (UniProt: `A0ABQ8MX81`, residues 333–372) emerged as the single highest-priority physiological candidate:

```
                          Serotransferrin Precursor (616 aa)
                                           │
                       ┌───────────────────┴───────────────────┐
                       ▼                                       ▼
             N-Lobe Transferrin Domain               C-Lobe Transferrin Domain
                       │
                       ▼ Residues 333–372
          [GVTWNNLKGR] KSCHTGLNRNAGWKVPDSAICGKTPDCTLYKFFSKGCAPG [ADPQSNMCEL]
          ──────────── ──────────────────────────────────────── ────────────
           N-Flank          Cryptic Bioactive AMP (40 aa)          C-Flank
         (Trypsin/Arg)                                          (Acidic/Endo)
```

#### Biochemical & Structural Attributes:
* **Sequence:** `KSCHTGLNRNAGWKVPDSAICGKTPDCTLYKFFSKGCAPG`
* **Net Charge:** $+4.0$ (5 Lysines [K1, K14, K23, K31, K35], 1 Arginine [R9], 1 Histidine [H4] vs. 2 Aspartates [D17, D26]).
* **Disulfide Architecture:** 4 Cysteine residues (`C3`, `C21`, `C27`, `C37`) predicted to form two internal disulfide staples that stabilize a compact cationic $\beta$-hairpin/loop architecture.
* **Amphiphilic Index:** $0.86$, conferring facial amphipathicity optimized for bacterial membrane insertion without destabilizing neutral zwitterionic eukaryotic erythrocyte membranes.
* **Aqueous Solubility:** High predicted solubility facilitating synthetic peptide handling and formulation.

#### Homology with Validated Ovotransferrin AMPs:
Sequence alignment against known antimicrobial peptide databases revealed remarkable structural conservation with validated ovotransferrin-derived peptides:

```text
L. rohita Cryptic AMP :  1  KSCHTGLNRNAGWKVPDSAICGKTPDCTLYKFFSKGCAPG  40
                            |||||||
P. sinensis sOT2 (Turtle):  1  KSCHTGL---------------------------------  ..  (Active AMP)
C. siamensis cOT2 (Croc) :  1  KSCHTGL---------------------------------  ..  (Active AMP)
```

The N-terminal heptapeptide `KSCHTGL` is identical to the active pharmacophore in **sOT2** (*Pelodiscus sinensis*) and **cOT2** (*Crocodylus siamensis*), which exhibit broad-spectrum bactericidal activity against both Gram-negative (*Aeromonas hydrophila*, *E. coli*) and Gram-positive (*Staphylococcus aureus*) pathogens.

---

### 4.3 Critical Appraisal: Discriminating True AMPs from Classifier Artifacts

A major insight from this proteome-scale screen is the necessity of **biological contextualization**:

```
                  Machine Learning Predictions (n = 4 Clusters)
                                        │
           ┌────────────────────────────┴────────────────────────────┐
           ▼                                                         ▼
Physiologically Authentic Candidate                       Classifier Artifacts / False Positives
  • Serotransferrin (Cluster 2)                            • GATA-2 Zinc Finger (Cluster 7)
    - Secreted host-defense protein                          - Nuclear transcription factor
    - Natural inflammatory cleavage                          - DNA-binding domain (C-X2-C...C-X2-C)
    - Conserved KSCHTGL active motif                       • TRAF3IP2 (Cluster 4) / CD1a (Cluster 6)
    - Proven antimicrobial homolog family                    - Intracellular / non-secreted
```

* **Cluster 7 (GATA-binding factor 2):** Sequence `...CVNCGATSTPLWRRDGTGHYLCNAC...` represents a classic canonical $\text{C}_2\text{H}_2/\text{C}_4$ zinc-finger DNA-binding motif. The spatial density of cysteine residues coupled with basic arginine/lysine residues artificially triggers CDP classification rules. Biologically, this nuclear transcription factor is never secreted or exposed to extracellular pathogens.
* **Cluster 4 (TRAF3IP2) & Cluster 6 (CD1a):** Represent intracellular signaling adaptors and membrane glycoprotein transmembrane segments that mimic amphipathic metrics in isolated fragments.
* **Conclusion:** Machine learning classifiers must be coupled with source-protein localization and physiological cleavage audits to prevent investment in false-positive artifacts.

---

## 5. Repository Organization

```text
labeo-rohita-amps/
├── Labeo_rohita_immune.fasta                        # Curated immune proteome (606 proteins, UniProt BAU-BD-2019)
├── file1.py                                         # Sliding window fragmentation script (15–50 aa)
├── file2.1.py                                       # Batch processor & Macrel prediction parser (CLP extraction)
├── file2.2.py                                       # Dynamic multi-family parser & threshold filter (Prob ≥ 0.90)
├── dedup_candidates.py                              # Union-find clustering & reverse source-protein tracing
├── final_results/                                   # Filtered, clustered candidate outputs
│   ├── combined_safe_amps_all_families.csv         # Consolidated high-confidence non-hemolytic hits
│   ├── combined_cdp_amps.csv                        # Filtered Cationic Disulfide-containing Peptides (CDP)
│   └── distinct_candidates_final.csv                # Final 4 non-overlapping candidate regions with source mapping
├── macrel_results_1/                                # Screening predictions for Partition 1
│   ├── README.md                                    # Macrel run metadata & citations
│   └── macrel.out.prediction.gz                     # Gzip-compressed Macrel predictions (~84k hits)
├── macrel_results_2/                                # Screening predictions for Partition 2
│   ├── README.md
│   └── macrel.out.prediction.gz                     # Gzip-compressed Macrel predictions (~83.6k hits)
├── macrel_results_3/                                # Screening predictions for Partition 3
│   ├── README.md
│   └── macrel.out.prediction.gz                     # Gzip-compressed Macrel predictions (~83.5k hits)
├── readme.txt                                       # Execution log & command history
├── .gitignore                                       # Git ignore rules for >100MB databases & caches
└── README.md                                        # Comprehensive project documentation
```

---

## 6. Execution, Reproducibility & Environment Setup

### 6.1 Computational Prerequisites
* Python 3.8+
* Conda / Mamba
* Bioinformatic Tools: `macrel`, `seqkit`, `biopython`, `pandas`

### 6.2 Step-by-Step Pipeline Execution

#### 1. Clone Repository & Set Up Environment
```bash
git clone https://github.com/Soash/labeo-rohita-amps.git
cd labeo-rohita-amps

# Create Conda environment with Macrel & SeqKit
conda create -y -n rohu_amp -c bioconda -c conda-forge macrel seqkit biopython pandas
conda activate rohu_amp
```

#### 2. Generate Sliding Window Fragment Library
```bash
# Generates Labeo_rohita_AMPs.fasta (~10M fragments, lengths 15–50)
python file1.py
```

#### 3. Partition Library for Parallel Screening
```bash
# Split the 10M fragments into 3 balanced parts
seqkit split2 -p 3 -O splits Labeo_rohita_AMPs.fasta
```

#### 4. Run Macrel Random Forest Screening
```bash
# Execute Macrel on each partition (adjust threads according to hardware)
macrel peptides --fasta splits/Labeo_rohita_AMPs.part_001.fasta --output ./macrel_results_1 --threads 28
macrel peptides --fasta splits/Labeo_rohita_AMPs.part_002.fasta --output ./macrel_results_2 --threads 28
macrel peptides --fasta splits/Labeo_rohita_AMPs.part_003.fasta --output ./macrel_results_3 --threads 28
```

#### 5. Filter & Consolidate Safe AMP Candidates
```bash
# Parses all macrel_results_*/ predictions, filters AMP_prob >= 0.90 & NonHemo
python file2.2.py
```

#### 6. Cluster Overlapping Windows & Trace Source Proteins
```bash
# Collapses overlapping window fragments via union-find & maps to source proteins
python dedup_candidates.py
```

---

## 7. Authorship & Correspondence

If you use this dataset, pipeline, or the identified candidate sequences, please cite:

> **Md Abu Sadat Soash** (2026). *In Silico Mining of Antimicrobial Peptide Candidates from the Immune Proteome of Labeo rohita Reveals a Cryptic Transferrin-Derived Peptide.* Department of Genetic Engineering and Biotechnology, Shahjalal University of Science and Technology, Sylhet, Bangladesh.  
> GitHub: [https://github.com/Soash/labeo-rohita-amps](https://github.com/Soash/labeo-rohita-amps)

### Inquiries
- **Md Abu Sadat Soash** — Department of Genetic Engineering and Biotechnology, Shahjalal University of Science and Technology, Sylhet, Bangladesh
- **Email:** [mdabu13@student.sust.edu](mailto:mdabu13@student.sust.edu)

---

## 8. References & Data Availability

1. **UniProt Reference Proteome:** *Labeo rohita* strain BAU-BD-2019 reference proteome available via [UniProt Proteomes (TaxID: 84645)](https://www.uniprot.org/taxonomy/84645).
2. **Macrel Prediction Engine:** Santos-Junior, C. D., et al. (2020). *Macrel: antimicrobial peptide screening in genomes and metagenomes.* **PeerJ**, 8, e10555. DOI: [10.7717/peerj.10555](https://doi.org/10.7717/peerj.10555).
3. **SeqKit Toolkit:** Shen, W., et al. (2016). *SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation.* **PLOS ONE**, 11(10), e0163962. DOI: [10.1371/journal.pone.0163962](https://doi.org/10.1371/journal.pone.0163962).
4. **Biopython:** Cock, P. J., et al. (2009). *Biopython: freely available Python tools for computational molecular biology and bioinformatics.* **Bioinformatics**, 25(11), 1422–1423. DOI: [10.1093/bioinformatics/btp163](https://doi.org/10.1093/bioinformatics/btp163).
5. **Ovotransferrin Antimicrobial Peptides:** Ibrahim, H. R., et al. (2000). *A structural phase of ovotransferrin for the emergence of antimicrobial activity.* **Biochimica et Biophysica Acta (BBA) - General Subjects**, 1523(2–3), 196–205.
