# Spatial-Proteomics-R-based-Data-Analysis-Pipeline

Analysis pipeline for DIA MS protein and peptide intensity datasets, which segregates/annotates protein categories of interest, carries out statistical analyses, and visualises results. 


<img src="https://github.com/user-attachments/assets/f73d7e31-4b80-426e-b872-4ed753d36b58" width="600"/>
<img src="https://github.com/user-attachments/assets/9ea0d4be-2f3b-45ba-b823-d9ee1898ce8e" width="600"/>
<img src="https://github.com/user-attachments/assets/44a9fea0-5ebe-447d-a5c4-aa8568680b3c" width="500"/>




## Data_filterer script
The first part of the pipeline (A) takes the protein and peptide datasets (either directly from DIA-NN or otherwise – any format works provided they have the necessary columns [details below]). The script accomplishes two main tasks: the formatting of peptide data and experiment feed ready for PLF analysis input, and the categorisation and filtering of datasets according to BM and ECM protein classes. This part essentially creates all datasets that may be wanted in later analysis, in one automatic step.

* Takes in protein and peptide datasets (any format works, providing requirements met: columns for UniProt AC, gene, and samples, plus a column for peptide sequence in peptide data).
* Input parameters (all input at the top of the page):
  * Experimental design:
    * Desired sample names (for graph annotations)
    * Experimental group names, with the number of samples per group
  * Select columns in datasets for: samples, gene, UniProt AC, peptide sequence (peptide data only)
  * Name folder for output datasets to be saved in (automatic)
  * Name an ‘experiment identifier’ – used to append output file names with, to prevent confusion when viewing data later
  * After this information is input, the code can be run in full, and all output files will be saved to a new folder (named in the parameters).
* The script first tidies the datasets by:
  * Removing rows with blank genes (non-human contaminants)
  * If the UniProt AC and/or gene columns contain multiple IDs, the first of these will be extracted and used for later steps (as per the UniProt website – the first AC is the primary, others are secondary)
* Creates labelled and filtered versions of the input protein and peptide data
  * Does this by searching protein/peptide UniProt ACs against the ECM and BM protein category lists
  * Annotates the original protein and peptide datasets with the categories, and creates datasets containing only ECM or only BM proteins (also with categories)
* Creates PLF-ready input peptide data and experiment feed, for all proteins, ECM proteins, and BM proteins
* Saves all datasets to a new folder


### ECM and BM protein lists:
MatrisomeDB, provided by The Matrisome Project, presents a comprehensive categorisation system for ECM proteins. The matrisomeDB may be searched here: https://matrisomedb.org/, but the full database can be downloaded as a spreadsheet containing 1027 genes (the homo sapiens master list) classified into matrisome divisions (core or associated) and six matrisome categories (glycoproteins, proteoglycans, collagens, regulators, secreted factors, and ECM-affiliated proteins). These matrisome categories are used to classify ECM proteins in this pipeline.

Similar to MatrisomeDB, BMbase is an online database of 222 human genes encoding BM zone proteins, divided into categories according to localisation, association, and experimental evidence. This database can be searched or downloaded – the latter was done for use in the pipeline. Only BMbase proteins with confirmed basement membrane zone localisation were used. The proteins were then subdivided into those with basement membrane subcellular localisation (termed ‘core BM proteins’), and the remainder (BM zone-associated proteins, such as cell membrane, secreted, or ECM proteins) (termed ‘confirmed BM proteins’). Therefore, the pipeline categorises BM proteins as core or confirmed, in addition to the MatrisomeDB categories (with overlap). 

Overall, any proteins with a BM category are segregated as BM proteins, likewise any protein with a matrisome category are segregated as ECM proteins, to produce separate datasets (as well as the full datasets with categories annotated).





## Ultimate_analysis script
The second part of the pipeline (B) performs all statistical analyses to produce final results and several data visualisations. Data undergoes filtering for peptides per protein and missing values, normalisation, imputation, then principal component analysis (PCA) and differential expression (DE) analysis on the protein data. Visualisations are generated for general data distribution, missing values, peptides and proteins per experimental group, and statistical analysis results.







