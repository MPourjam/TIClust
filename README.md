# TIC (Taxonomy Informed Clustering) Pipeline 

## Overview

The **Taxonomy Informed Clustering (TIC) pipeline** is a novel approach for processing 16S rRNA amplicon datasets in diversity analyses. TIC leverages classifier-assigned taxonomy to refine clustering, ensuring that only sequences sharing the same taxonomic path are grouped together [1]. This method is particularly useful for analyzing bacterial diversity using 16S rRNA gene amplicon sequencing data.

Inspired by the initial TIC pipeline as described in [Taxonomy Informed Clustering, an Optimized Method for Purer and More Informative Clusters in Diversity Analysis and Microbiome Profiling](https://doi.org/10.3389/fbinf.2022.864597), this project aims to provide a simple and easy-to-use version of the original TIC.


## Key Features

*   **Taxonomy-Driven Clustering**: TIC taxonomically classifies each sequence before clustering, using the taxonomic information to guide and constrain the clustering process. This approach divides the dataset into subsets based on assigned taxonomies, preventing the merging of sequences from diverse lineages.
*   **Modular Design**: The pipeline is designed with a modular structure, making it easy to modify and extend.
*   **Comprehensive Pipeline**: TIC offers a complete, automated pipeline for diversity analyses, from raw reads to compositional tables.

## Pipeline Steps

The TIC pipeline involves several key steps:

1.  **Filling upto Order Level**: All bacterial sequences with taxonomy not known upto order level are filled with `NA__<level>__<parent>`
2.  **Complete Family Levels**: Sequences of orders get clustered at *0.90* sequences similarity using *uclust* and assigned their family level taxonomy (e.g: *FOTU11*)
3.  **Complete Genus Levels**: Sequences of families get clustered at *0.95* sequence similarity using *uclust* and assigned their genus level taxonomy (e.g: *GOTU11*)
4.  **Complete Species Levels**: Sequences of genera get clustered at *0.987* similarity using *uclust* and assigned their family level taxonomy (e.g: *SOTU11*)


## Installation
To install the TIC pipeline, follow these steps:

### PIP
```bash
pip install tic
```

## Run TIC
> Currently simple version of TIC is available in this repo

After tic is installed you can run *simple tic* by commands below:

```bash
simpletic -f <path/to/taxed_fasta_file.fasta>
```

or 

```bash
python3 -m simpletic -f <path/to/taxed_fasta_file.fasta>
```

Sequences in the input fasta file should have their taxonomies as **last** part of their sequence header starting with **tax=**. For example:
`>Zotu1 some other header parts tax=Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia-Shigella`.
Note that taxonomies levels are separated by semicolon (`;`).


## Outputs

TIC will create a directory (i.e: `TIC-WD`) in the same directory of the input Fasta file. Inside the directory resides files listed below:

1. `TIC-FullTaxonomy.fasta`: *bacterial* sequences with full-level taxonomy filled by TIC
2. `Map-FOTU-GOTU.tab`: two columns 'FOTU' and 'GOTU' which maps genus-level clusters (i.e: GOTUs or gOTUs) to their parent family-level clusters (i.e: FOTUs or fOTUs)
3. `Map-GOTU-SOTU.tab`: two columns 'GOTU' and 'SOTU' which maps species-level clusters (i.e: SOTUs or sOTUs) to their parent genus-level clusters (i.e: GOTUs or gOTUs)
4. `Map-SOTU-ZOTU.tab`: two columns 'SOTU' and 'ZOTU' which maps zero-radius operational taxonomic units (i.e: ZOTUs or zOTUs) to their parent species-level clusters (i.e: SOTUs or sOTUs)
5. `Non-Bacteria-Sequences.fasta`: all non-bacterial sequences in input Fasta file.


## Advantages

TIC enhances the clustering process by utilizing taxonomic information acquired beforehand, leading to higher cluster quality and purity. This approach also enables the proper placement of unassigned sequences within the taxonomy.

## Availability

TIC as a python package (i.e: current repository) is available on [SimpleTIC](https://github.com/MPourjam/SimpleTIC).
The initial TIC pipeline as described in [Taxonomy Informed Clustering, an Optimized Method for Purer and More Informative Clusters in Diversity Analysis and Microbiome Profiling](https://doi.org/10.3389/fbinf.2022.864597) is freely available on [GitHub](https://github.com/Lagkouvardos/TIC) [1].
