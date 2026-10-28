conda activate base

pip install biopython 



conda create -y -n amplify python=3.6

conda activate amplify

conda install -y -c bioconda amplify

conda deactivate

AMPlify --help

AMPlify -m imbalanced -s Rohu_AMPlify_Library.fasta -od ./rohu_results -of tsv

AMPlify -m imbalanced -s Rohu_AMPlify_Library_test.fasta -od ./rohu_results -of tsv



seqkit split2 -p 3 -O splits Labeo_rohita_AMPs.fasta

conda create -n macrel -c bioconda -c conda-forge macrel -y

conda activate macrel

macrel peptides --fasta splits/Labeo_rohita_AMPs.part_001.fasta --output ./macrel_results_1 --threads 28

macrel peptides --fasta Labeo_rohita_AMPs.fasta --output ./macrel_results_1 --threads 28



Generated 10049499 unique fragments (lengths 15-50)!

(base) bionanosense@bionanosense:~/soash/labeo_rohita_amps$ seqkit split2 -p 3 -O splits Labeo_rohita_AMPs.fasta
[INFO] split seqs from Labeo_rohita_AMPs.fasta
[INFO] split into 3 parts
[INFO] write 3349833 sequences to file: splits/Labeo_rohita_AMPs.part_001.fasta
[INFO] write 3349833 sequences to file: splits/Labeo_rohita_AMPs.part_002.fasta
[INFO] write 3349833 sequences to file: splits/Labeo_rohita_AMPs.part_003.fasta
(base) bionanosense@bionanosense:~/soash/labeo_rohita_amps$ 







