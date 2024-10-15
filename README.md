# 🖥️🧬 Gene Prediction

## 🔍 Overview

Genes are subsequences within transcripts that can be translated into proteins by ribosomes. They contain a reading frame made up of consecutive triplets, starting with an initiation codon (e.g., 'AUG', 'UUG', 'CUG', 'AUU', or 'GUG') and ending with a stop codon (either 'UAA', 'UAG', or 'UGA'). All codons must be in the same reading frame.

Upstream of the initiation codon, there is typically a motif that aids in translation initiation by allowing the binding of the 16S ribosomal RNA subunit. This sequence, known as the Shine-Dalgarno sequence (AGGAGGUAA), helps position the ribosome correctly [\[Shine and Dalgarno, 1973\]](https://www.sciencedirect.com/science/article/pii/0022283673905287). It is worth noting that the Shine-Dalgarno motif doesn't necessarily have to be in the same reading frame as the initiation codon, and it may sometimes be incomplete.

![Gene Prediction](assets/gene-prediction.png)

While many organisms still lack experimentally verified genome annotations, gene prediction remains crucial for the automatic annotation of genomes. Several tools and methods are available for this task [\[List of gene prediction software\]](https://en.wikipedia.org/wiki/List_of_gene_prediction_software).

In this project, we aim to develop a simple approach to predict prokaryotic genes by detecting reading frames and identifying the Shine-Dalgarno motif. Our objective is to predict genes for the reference genome of *Listeria monocytogenes* EGD-e (sequenced by the Institut Pasteur), which has 2,867 known genes. You can find the genome [here](https://www.ncbi.nlm.nih.gov/genome/browse/#!/proteins/159/159660%7CListeria%20monocytogenes%20EGD-e/).

## 🔄Installation

To set up the environment and install the required dependencies, use the following commands:

```bash
conda env create -f environment.yml
conda activate gene-prediction
```

## 🧑‍💻️Usage

First, clone the repository and navigate to the project folder:

```bash
git clone git@github.com:zhukovanadezhda/gene-prediction.git
cd gene-prediction
```

To run the gene prediction script, use the following command:

```
python3 gpred/gpred.py -i data/listeria.fna \
                       -p results/predicted_gene_positions.csv \
                       -o results/predicted_genes.fasta
```

The following command will analyze the Listeria genome file located at `data/listeria.fna` to predict genes. It specifies a minimum gene length of 50 nucleotides, a maximum allowable distance of 16 nucleotides between the start codon and the Shine-Dalgarno sequence, and requires at least a 40-nucleotide gap between consecutive genes. The predicted gene sequences and their corresponding positions will be saved in `results/predicted_genes_positions.csv` and `results/predicted_genes.fasta`, respectively.

### ⚙️ Command-line Options:

```
  -h, --help                      Show help message and exit
  -i GENOME_FILE                  Complete genome file in FASTA format
  -g MIN_GENE_LEN                 Minimum gene length to consider (default: 50)
  -s MAX_SHINE_DALGARNO_DISTANCE  Max distance from start codon to search for Shine-Dalgarno motif (default: 16)
  -d MIN_GAP                      Minimum gap between two genes (excluding Shine-Dalgarno box, default: 40)
  -p PREDICTED_GENES_FILE         Output CSV file with predicted gene positions
  -o FASTA_FILE                   Output FASTA file containing predicted gene sequences
```

## 🎁 Test the results 

To assess the accuracy of our gene predictions, we can compare them against a reference set of known genes using [jvenn](https://jvenn.toulouse.inra.fr/app/example.html), an online tool developed by INRA for Venn diagram analysis. The comparison involves the following datasets:

1. **Predicted genes**: Genes identified by our program.
2. **Reference genes (Prodigal)**: From `data/prodigal.csv`, which contains genes predicted by the Prodigal software.
3. **Positions of reference genes**: From `data/positions.csv`, detailing the location of known genes within the Listeria genome.

![Comparison](assets/venn_chart.png)  

This comparison reveals that relying solely on the Shine-Dalgarno motif for gene prediction is insufficient to capture all genes accurately. However, it can still help identify certain genes that Prodigal does not predict.


