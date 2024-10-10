# Gene Prediction

Genes correspond to subsequences of transcripts that can be translated into proteins by the ribosome. They contain a reading frame made up of consecutive triplets starting from an initiation codon ('AUG', 'UUG', 'CUG', 'AUU', or 'GUG') and ending with a stop codon ('UAA', 'UAG', or 'UGA'). These codons are within the same reading frame!

![Gene Prediction](assets/gene-prediction.png)

Upstream of the initiation codon, there's a motif that enables the initiation of translation by binding the 16S subunit of the ribosomal RNA: the **Shine-Dalgarno sequence** (AGGAGGUAA) [Shine and Dalgarno 1973]. This motif doesn't have to be in the same reading frame as the initiation codon and can be incomplete.

Currently, only a few organisms have experimentally verified gene annotations. Therefore, gene prediction remains a crucial task for the automatic annotation of genomes. Many software tools and approaches exist for this task:
[Gene Prediction Software](https://en.wikipedia.org/wiki/List_of_gene_prediction_software)

In this project, we develop a simple approach to predict genes in prokaryotes, based on detecting reading frames and the Shine-Dalgarno motif. The goal of this project is to predict genes in the reference genome of **Listeria monocytogenes EGD-e** (assembled and sequenced by the Pasteur Institute), which contains 2867 genes:
[NCBI Listeria monocytogenes EGD-e Genome](https://www.ncbi.nlm.nih.gov/genome/browse/#!/proteins/159/159660%7CListeria%20monocytogenes%20EGD-e/)
