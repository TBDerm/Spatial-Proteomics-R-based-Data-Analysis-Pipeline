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


* Takes in protein and peptide datasets (including annotated protein categories from the previous step of the pipeline, if desired)
* Input parameters (all input at the top of the page):
  * Experimental design:
     * Desired sample names (for graph annotations)
    * Experimental group names, with the number of samples per group
  * Select columns in datasets for: samples, gene, UniProt AC, peptide sequence (peptide data only)
  * Name folder for output datasets to be saved in (automatic)
  * Name an ‘experiment identifier’ – used to append output file names with, to prevent confusion when viewing data later
  * Select type of normalisation wanted (options justified in section 3.1.4)
    * No normalisation
    * Median absolute deviation normalisation
    * Variance stabilising normalisation
    * Median centring normalisation
    * Cyclic loess normalisation
    * Robust linear scaling
    * Median-based scaling
    * Median centring the median absolute deviation normalisation
    * Median centring then standard deviation (or z-score) normalisation
  * Select type of imputation wanted (options justified in section 3.1.4)
    * Predictive mean matching
    * Classification and regression trees
    * Lasso linear regression
    * K nearest neighbours (K=3)
  * Choose colours for the colour scheme of all output visualisations (one colour per experimental group, and high-low colours for heatmaps)
* After this information is input, the code can be run in full, and all output files will be saved to a new folder (named in the parameters).
* Generates venn diagrams of proteins, and peptides, in each experimental group (fig. 7a)
  * This shows the number of proteins/peptides identified in each experimental condition, and overlaps between conditions 
  * Also generates a table of proteins within each intersection of the venn 
* Using the peptide dataset, counts the number of unique peptides per protein, and amends this to the protein dataset as a new column
* Filters protein dataset to exclude any proteins with <2 peptides
  * Generates a missing values plot (fig. 7b) showing proportion of missing values per sample, prior to filtering for missing values
* Filters protein dataset to exclude any proteins with missing values in >40% of samples in any experimental group (therefore all remaining proteins have  values present in >60% of samples in every experimental groups)
  * Generates a second missing values plot showing proportion of missing values per sample, after filtering for missing values
* Remaining proteins are termed ‘valid’ (proteins with >1 peptide and >60% present values per experimental group) – peptide data is filtered to include only peptides for valid proteins
  * The filtered ‘valid’ protein and peptide datasets are exported
* Log2 -transforms the protein sample data
* Normalises the data, using user-selected method
  * The normalised protein data is exported
* Imputes the data, independently per experimental group, using user-selected method
  * The imputed protein data is exported
* Generates a variety of data visualisations (selection shown in figure 7) (based on the filtered imputed protein data)
  * Heatmap of all proteins in all samples, divided into experimental groups
  * Heatmap of the top 30 proteins by median intensity value
  * Z-scored heatmap of the top 30 proteins by median intensity value (fig. 7d)
  * Violin and boxplots to show protein data distribution per sample, divided into experimental groups (fig. 7h)
  * Frequency distribution curves of protein intensity values, per experimental group and per sample (fig. 7f) – these allow the user to check that data distribution is normal and similar for each sample/group after normalisation
  * Frequency distribution histograms of protein intensity values, per experimental group (fig. 7g) – an alternative way of displaying the above
  * Histograms of number of peptides per protein, across all samples and per experimental group (fig. 7c) – informs the user how many peptides were identified for the majority, or minority, of proteins
* Performs PCA on the imputed protein data
  * Generates a large PCA plot of PC1 and PC2, and a composite of all PC axes combinations up to PC4 (fig. 7e) – shows any separation of experimental groups, thus attributing a percentage of variance in the data (derived from PCs, hence the different axes display) to differences between these experimental groups
  * Built-in options to customise PCA plots: with/without frames, with/without 95% confidence ellipses (default is with frames, no ellipses)
* Performs DE analysis on the imputed protein data
  * Uses the limma R package system, which fits a linear model and empirical bayes moderation to generate pairwise comparisons. P-values are corrected for multiple comparisons using Benjamini Hochberg.
  * DE results (log2-fold change (FC) values, t and B values, P and P-adjusted (P-adj) values) are appended to the full protein dataset (with categories if imported) and full results data exported
  * Generates volcano plots (fig. 7i) (one per combination pair of experimental groups, if there are over two), with and without protein gene labels for proteins beyond the significance threshold
  * Built-in options to customise volcano plots: p-adj value cutoff, fold change cut off, and colour of datapoints beyond these cutoffs (considered significant). Default settings are P-adj<0.05, FC>1.5, and significant points coloured red
* All datasets and results tables are saved as csv files in a new folder. All visualisations are saved in a PDF file (one plot per page)






