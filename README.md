# Shotgun-Metagenomics-Assignment3
Shotgun metagenomic analysis comparing gut microbiomes of vegans and omnivores

# Introduction

The human gut microbiome is a complex microbial community that plays an essential role in human health. Gut microorganisms perform important metabolic functions, including the biosynthesis of vitamins and the breakdown of otherwise indigestible compounds, while also interacting with the host through the production of metabolites that can have beneficial or detrimental effects [1]. Disruptions in the balance of these microbial communities, commonly referred to as dysbiosis, have been associated with numerous diseases such as obesity, inflammatory bowel disease, cardiovascular disease, neurological disorders, and certain cancers [1].

Diet is one of the most important factors influencing the composition and function of the gut microbiome. Diets rich in plant-based foods provide high levels of fiber and complex carbohydrates that gut microbes ferment to produce short-chain fatty acids (SCFAs), which help maintain gut barrier integrity and support immune regulation [2]. Plant-based foods also contain polyphenols that promote beneficial bacterial populations and help prevent inflammation and pathogen colonization [2]. In contrast, diets high in animal-derived foods increase protein fermentation in the gut, which may lead to inflammation, reduced SCFA production, and the formation of metabolites such as trimethylamine N-oxide (TMAO), which has been associated with cardiovascular disease and colorectal cancer [2].

Dietary patterns have also been linked to shifts in specific microbial taxa. For example, the genus Prevotella is commonly associated with plant-rich, agrarian diets, while Bacteroides species are more frequently linked to diets high in animal protein and fat [1]. Changes in dietary patterns, particularly the adoption of Westernized diets characterized by high fat and protein intake and reduced fiber consumption, have been associated with reduced microbial diversity and potential negative health outcomes [1].

Advances in metagenomic sequencing have greatly expanded our understanding of the gut microbiome and its relationship with diet and health. Shotgun metagenomics allows researchers to analyze all genetic material within a sample, enabling species- and strain-level characterization of microbial communities and providing insights into their functional potential [1]. Because individual microbial species and strains can possess distinct genomic and metabolic capabilities, this approach provides a more comprehensive understanding of microbiome structure and function than traditional taxonomic profiling alone [1].

Shotgun metagenomic sequencing provides a comprehensive approach for studying microbial communities. Although it is more expensive than targeted approaches such as 16S rRNA gene sequencing, it offers several important advantages [3]. In addition to identifying the taxonomic composition of microbial communities, shotgun metagenomics also enables the analysis of genes present in the microbiome, allowing researchers to investigate the functional capabilities of microbial populations [3]. Numerous bioinformatic tools have been developed to analyze shotgun sequencing data, including methods for taxonomic classification, prediction of gene abundance such as antibiotic resistance or virulence genes, reconstruction of metagenome-assembled genomes (MAGs), and metabolic pathway analysis [3]. Many of these tools can directly access sequencing data deposited in public repositories using accession numbers, which simplifies large-scale analysis and facilitates reproducible microbiome research.

Assessing sequencing quality is an essential first step in metagenomic analysis, as downstream taxonomic classification and functional profiling depend on high-quality input data. Poor-quality reads can lead to inaccurate assignments and biased results, making quality control a critical component of the pipeline. To evaluate read quality, the tool FastQC will be used, which provides comprehensive visual summaries of sequencing quality metrics, including per-base quality scores, sequence composition, and potential technical biases [4]. FastQC also identifies the encoding of base quality scores in FASTQ files, ensuring correct interpretation of quality values [4]. Overall, this step ensures that the data are suitable for reliable downstream analysis.

FASTQ preprocessing will be performed using fastp, a fast and comprehensive tool for quality control and read filtering [5]. High-quality preprocessing is essential in metagenomic analyses, as low-quality reads and adapter contamination can negatively impact downstream taxonomic classification. Compared to other commonly used tools such as Trimmomatic and Cutadapt, fastp provides an integrated solution that performs multiple functions within a single step, including adapter trimming, quality filtering, length filtering, and generation of quality control reports [5]. In contrast, Trimmomatic and Cutadapt typically require additional tools such as FastQC or MultiQC to assess data quality before and after processing, increasing computational time and resource usage [5]. Fastp streamlines this workflow by combining preprocessing and quality assessment while maintaining high efficiency and speed [5]. Its ability to handle multiple FASTQ files in parallel and generate detailed HTML reports further enhances its usability [5]. Therefore, fastp was selected for its efficiency, comprehensive functionality, and ability to simplify the preprocessing pipeline without compromising performance. 

Taxonomic classification of metagenomic reads will be performed using Kraken2, a fast and memory-efficient tool that assigns reads based on k-mer matching to reference genomes [6]. Compared to its predecessor Kraken1, Kraken2 significantly reduces memory requirements while maintaining high accuracy and increasing classification speed [6]. This improvement allows the use of larger and more comprehensive reference databases, which enhances taxonomic resolution without requiring extensive computational resources [6]. Additionally, Kraken2 remains compatible with downstream tools such as Bracken for species-level abundance estimation [6]. Due to its improved efficiency, scalability, and accuracy over Kraken1, Kraken2 is well suited for large-scale metagenomic analyses. A confidence threshold of 0.1 will applied in Kraken2 to balance classification accuracy and read retention. Increasing the confidence threshold reduces false positive assignments but can result in a higher proportion of unclassified reads [7]. As this study involves complex microbial communities with potentially high diversity, a moderate threshold was selected to maintain sensitivity while still limiting spurious classifications [7]. A reduced Kraken2 reference database (~16 GB) will used due to computational limitations associated with the full database (>100 GB). While larger databases and higher confidence thresholds can improve classification precision, practical considerations such as computational resources and the need to retain sufficient read coverage require a balanced approach.

Following taxonomic classification with Kraken2, Bracken will be used to estimate taxonomic abundances. Bracken (Bayesian Reestimation of Abundance after Classification with KrakEN) refines Kraken2 output by re-estimating read counts at a specified taxonomic level, such as species or genus, using information about genome length and k-mer distributions [8]. This approach improves abundance estimation accuracy, particularly in cases where closely related species share similar sequences and may be misclassified by read-level methods alone [8]. By redistributing reads assigned at higher taxonomic levels, Bracken provides more reliable species-level abundance estimates [8]. Additionally, Bracken applies a minimum read threshold to reduce false positives, ensuring that only taxa with sufficient support are retained [8]. Therefore, Bracken will be used to enhance the accuracy and reliability of taxonomic abundance estimates in this study.

The phyloseq package will be used to facilitate microbiome data analysis by integrating abundance data, taxonomic information, and sample metadata within a single framework [9]. This structured approach simplifies data management and ensures consistency across analyses, as all components remain linked throughout processing [9]. Phyloseq also provides convenient functions for data transformation, visualization, and statistical analysis, enabling efficient and reproducible workflows [9]. Its ability to organize and handle complex sequencing data makes it well suited for comparative microbiome studies.

For differential abundance analysis, ANCOM-BC2 will be used due to its improved performance over the original ANCOM-BC method [10]. Previous studies have shown that ANCOM-BC2 provides better control of the false discovery rate (FDR) while maintaining high statistical power, particularly in complex or edge-case scenarios where earlier methods may perform less reliably [10]. These improvements make ANCOM-BC2 more robust for identifying differentially abundant taxa in microbiome datasets [10]. Additionally, ANCOM-BC2 is well suited for small to moderate sample sizes, making it appropriate for this study [10]. Therefore, ANCOM-BC2 was selected to ensure accurate and reliable detection of differential abundance between groups..

In this study, shotgun metagenomic sequencing data from human gut samples will be analyzed to compare microbial communities between individuals following vegan and omnivorous diets. This analysis aims to investigate how dietary patterns influence gut microbiome composition and to explore potential differences in microbial taxa associated with plant-based and animal-based diets.

# Methods

## Computational Environment

All analyses were performed on a high-performance computing (HPC) environment using Linux command-line tools. Bash scripts were used to automate data download, preprocessing, and taxonomic classification steps to ensure reproducibility. Downstream statistical analyses were conducted in R (v4.5.1) using relevant Bioconductor and CRAN packages.

## Data Acquisition and Quality Control

Shotgun metagenomic sequencing data were obtained from the NCBI Sequence Read Archive (SRA) using accession numbers SRR8146935 (omnivore), SRR8146936 (omnivore), SRR8146938 (omnivore), SRR8146951 (vegan), SRR8146952 (vegan), and SRR8146954 (vegan). Data were downloaded using the SRA Toolkit (v3.0.9) with the `prefetch` command, followed by conversion to FASTQ format using `fasterq-dump` [11,12]. The resulting paired-end FASTQ files were compressed using `gzip` to reduce storage requirements and improve computational efficiency.

Initial quality assessment of sequencing reads was performed using FastQC (v0.12.1) [13]. Each FASTQ file was analyzed individually to evaluate per-base sequence quality, GC content, sequence duplication levels, and adapter contamination.

## Read Preprocessing

Quality control and preprocessing were conducted using fastp (v1.0.1) [14]. Paired-end reads were processed using default parameters with automatic adapter detection enabled. Low-quality bases were trimmed, and reads with an average Phred score below 20 or length shorter than 50 bp were removed. Fastp was executed with multi-threading enabled (`-w`) to improve computational efficiency. HTML and JSON reports were generated for each sample to summarize preprocessing statistics.

## Taxonomic Classification

Taxonomic classification was performed using Kraken2 (v2.1.6) [15]. A pre-built standard reference database (~16 GB) was downloaded and used for classification [16]. Reads were classified using paired-end mode with a confidence threshold of 0.1 (`--confidence 0.1`) to balance classification accuracy and read retention. Kraken2 assigns taxonomic labels based on exact k-mer matches to reference genomes. Output included classification files and summary reports (`--report`) containing taxonomic abundance estimates.

## Abundance Estimation

Species-level abundance estimation was refined using Bracken (v3.0) [17]. Bracken was run using the same Kraken2 database with a read length parameter of 150 bp (`-r 150`) and species-level classification (`-l S`). Kraken2 report files were used as input, and Bracken re-estimated read counts to correct for biases introduced by shared k-mers among closely related taxa.

## Data Import and Processing in R

Bracken output files were imported into R and combined into a single dataset. Species-level abundance tables were constructed using estimated read counts (`new_est_reads`) and converted into a sample-by-species matrix. Data were analyzed using the phyloseq package (v1.52.0), which integrates abundance data, taxonomy, and metadata within a unified object structure [18].

## Rarefaction Analysis

Rarefaction curves were generated using the `rarecurve()` function from the vegan package (v2.7.3) to assess sequencing depth and sampling completeness [19]. Raw count data were used, and the abundance matrix was transposed to ensure samples were represented as rows. Curves were generated with a step size of 5000 reads.

## Relative Abundance Analysis

Raw counts were converted to relative abundances using `transform_sample_counts()` in phyloseq [18]. Data were reshaped using `psmelt()` and visualized using ggplot2 (v4.0.2) [18,20]. The top 10 most abundant species were identified, and remaining taxa were grouped into an “Other” category. Relative abundance was visualized using stacked bar plots, faceted by diet group.

## Alpha Diversity Analysis

Alpha diversity was calculated using the Shannon index via `estimate_richness()` in phyloseq [18]. Samples were grouped by diet (omnivore vs. vegan), and diversity distributions were visualized using boxplots with overlaid jittered points. Statistical differences between groups were assessed using a Welch two-sample t-test.

## Beta Diversity Analysis

Beta diversity was assessed using Bray–Curtis dissimilarity calculated with `distance()` in phyloseq [18]. Principal Coordinates Analysis (PCoA) was performed using `ordinate(method = "PCoA")`[18]. Ordination results were visualized using ggplot2, with samples colored by diet group. Statistical significance of group differences was evaluated using PERMANOVA via `adonis2()` in the vegan package [19].

## Differential Abundance Analysis

Differential abundance analysis was conducted using ANCOMBC (ANCOM-BC2) (v2.13.1) [21]. Species-level abundances were analyzed using diet as the main variable (`fix_formula = "Diet"`). P-values were adjusted for multiple testing using the Holm method, and taxa with adjusted q-values < 0.05 were considered statistically significant. Log fold changes and standard errors were extracted and visualized using ggplot2.

## Statistical Analysis and Visualization

All statistical analyses and visualizations were performed in R using phyloseq, vegan, ANCOMBC, and ggplot2. Figures were exported using `ggsave()` or base R graphics functions for inclusion in the final report.


# Results
<img width="3000" height="1800" alt="image" src="https://github.com/user-attachments/assets/bdddff8a-7347-4afe-8b9d-a9d3e50fb1fe" />

**Figure 1. Relative abundance of microbial species across omnivore and vegan samples.**
Relative abundance of the top 10 most abundant species is shown for each sample, with remaining taxa grouped as “Other.” Samples are grouped by diet (omnivore vs. vegan) and displayed as stacked bar plots. Each bar represents one sample, and colors indicate individual species contributions to the overall community composition.

Relative abundance analysis was performed to assess microbial community composition across samples (Figure 1). Across both diet groups, a small number of dominant species accounted for the majority of the community, with remaining taxa grouped as “Other.” Omnivore samples exhibited greater variability in species composition, with some samples showing higher contributions from Faecalibacterium prausnitzii and others dominated by different taxa. In contrast, vegan samples showed a more consistent pattern, with Prevotella copri representing a substantial proportion of the microbial community in multiple samples. Despite these trends, considerable variability was observed between individual samples within each diet group, indicating that inter-individual differences may contribute substantially to overall microbial composition.

<img width="2100" height="2100" alt="image" src="https://github.com/user-attachments/assets/c3a67058-97be-4ce8-8a1c-e846ca185a05" />

**Figure 2. Alpha diversity (Shannon index) of microbial communities across diet groups.**
Alpha diversity was measured using the Shannon index for each sample. Samples are grouped by diet (omnivore vs. vegan). Boxplots display the distribution of diversity values within each group, with individual samples shown as overlaid points.

Alpha diversity was assessed using the Shannon index to evaluate within-sample microbial diversity (Figure 2). Omnivore samples exhibited a slightly higher median Shannon diversity compared to vegan samples, although variability was observed within both groups. Vegan samples showed a narrower distribution of diversity values, indicating more consistent diversity across samples. 

**Table 1. Comparison of alpha diversity (Shannon index) between diet groups.**
Mean Shannon diversity values are shown for omnivore and vegan samples. A Welch two-sample t-test was performed to assess differences in alpha diversity between groups.
| Group | Mean Shannon | SD | n |
| :---------: | :-----------: | :----: | :---: |
| Omnivore | 2.869663 | 0.6500763 | 3 |
| Vegan | 2/673013   | 0.2538101 | 3 |

Welch two-sample t-test: t = 0.49, df = 2.60, p = 0.6637
95% confidence interval: −1.21 to 1.60

Alpha diversity was compared between diet groups using a Welch two-sample t-test (Table 1). The mean Shannon diversity was slightly higher in omnivore samples (mean = 2.87) compared to vegan samples (mean = 2.67). However, this difference was not statistically significant (t = 0.49, df = 2.60, p = 0.6637), and the 95% confidence interval (−1.21 to 1.60) included zero, indicating no evidence of a difference in within-sample diversity between diet groups. The Welch two-sample t-test revealed no significant difference in alpha diversity between diet groups, suggesting that within-sample microbial diversity is comparable between omnivore and vegan individuals in this dataset. 

<img width="2100" height="2100" alt="image" src="https://github.com/user-attachments/assets/946d2385-3026-4cf1-a992-2157e147980f" />

**Figure 3. Principal Coordinates Analysis (PCoA) of microbial community composition based on Bray–Curtis dissimilarity.**
PCoA was performed using Bray–Curtis distances calculated from species-level relative abundance data. Samples are colored by diet group (omnivore vs vegan), and each point represents one sample.

Beta diversity was assessed using Bray–Curtis dissimilarity and visualized using Principal Coordinates Analysis (PCoA) to evaluate differences in microbial community composition between samples (Figure 3). Samples showed partial clustering by diet group, with omnivore samples generally occupying a distinct region of the ordination space compared to vegan samples. However, some overlap between groups was observed, indicating that differences in community composition are not completely separated by diet. Additionally, variation within diet groups was evident, with one vegan sample (SRR8146954) appearing more distant from other samples along the second principal coordinate. Overall, these results suggest moderate differences in microbial community structure between diet groups, with substantial inter-individual variability.

**Table 2. PERMANOVA results testing differences in microbial community composition between diet groups.**
Permutational multivariate analysis of variance (PERMANOVA) was performed using Bray–Curtis dissimilarity to assess differences in community composition between omnivore and vegan samples.
| Source |	Df | Sum of Squares	| R² | F	| p-value |
| :---------: | :-----------: | :----: | :---: | :---: | :---: |
| Diet |	1	| 0.23389 |	0.20835 |	1.0528 |	0.20 |
| Residual | 4 |	0.88866 |	0.79165 |	— | — |
| Total |	5	| 1.1255 |	1.00000 |	— | — |

Permutational multivariate analysis of variance (PERMANOVA) was performed to assess whether microbial community composition differed between diet groups based on Bray–Curtis dissimilarity (Table 2). Diet explained approximately 20.8% of the variation in community composition (R² = 0.21); however, this effect was not statistically significant (F = 1.05, p = 0.20). These results indicate that, although some variation in microbial composition is associated with diet, the differences between omnivore and vegan groups are not strong enough to be distinguished from within-group variability in this dataset.

<img width="2100" height="2100" alt="image" src="https://github.com/user-attachments/assets/bcb6e361-a063-45da-8440-1f7d941d70b6" />

**Figure 4. Differential abundance analysis of microbial species between diet groups using ANCOM-BC2.**
Log fold changes (vegan vs omnivore) are shown for the top 15 taxa ranked by absolute effect size. Points represent estimated log fold changes, and horizontal bars indicate standard errors. The red vertical line at zero indicates no difference between groups.

Differential abundance analysis was performed using ANCOM-BC2 to identify taxa associated with diet (Figure 4). No taxa were found to be significantly different between omnivore and vegan groups after multiple testing correction (q < 0.05). Although several taxa exhibited variation in log fold change, with some enriched in vegan samples (positive values) and others enriched in omnivore samples (negative values), all confidence intervals overlapped zero. These results indicate that no individual species showed a statistically significant association with diet in this dataset, despite observable trends in relative abundance.

<img width="1392" height="860" alt="image" src="https://github.com/user-attachments/assets/94a4999b-02e8-4e9f-858b-6b81fca7bb64" />

**Figure 5. Rarefaction curves of microbial species richness across samples.**
Rarefaction curves were generated using species-level count data to assess sequencing depth and sampling completeness. Each curve represents one sample, showing the relationship between sequencing depth (number of reads) and observed species richness.

Rarefaction analysis was performed to evaluate whether sequencing depth was sufficient to capture microbial diversity (Figure 5). All samples exhibited curves that approached a plateau, indicating that additional sequencing would likely yield minimal increases in observed species richness. While differences in total richness were observed between samples, the overall shape of the curves suggests adequate sampling depth across all samples.

# Discussion


# References

[1] De Filippis, F., Pasolli, E., Tett, A., Tarallo, S., Naccarati, A., De Angelis, M., Neviani, E., Cocolin, L., Gobbetti, M., Segata, N., & Ercolini, D. (2019). Distinct Genetic and Functional Traits of Human Intestinal Prevotella copri Strains Are Associated with Different Habitual Diets. Cell Host & Microbe, 25(3), 444-453.e3. https://doi.org/10.1016/j.chom.2019.01.004

[2] Fackelmann, G., Manghi, P., Carlino, N., Heidrich, V., Piccinno, G., Ricci, L., Piperni, E., Arrè, A., Bakker, E., Creedon, A. C., Francis, L., Capdevila Pujol, J., Davies, R., Wolf, J., Bermingham, K. M., Berry, S. E., Spector, T. D., Asnicar, F., & Segata, N. (2025). Gut microbiome signatures of vegan, vegetarian and omnivore diets and associated health outcomes across 21,561 individuals. Nature Microbiology, 10(1), 41–52. https://doi.org/10.1038/s41564-024-01870-z

[3] Yan, J., Liao, C., Taylor, B. P., Fontana, E., Amoretti, L. A., Wright, R. J., Littmann, E. R., Dai, A., Waters, N., Peled, J. U., Taur, Y., Perales, M.-A., Siranosian, B. A., Bhatt, A. S., van den Brink, M. R. M., Pamer, E. G., Schluter, J., & Xavier, J. B. (2022). A compilation of fecal microbiome shotgun metagenomics from hematopoietic cell transplantation patients. Scientific Data, 9(1), Article 219. https://doi.org/10.1038/s41597-022-01302-9

[4] Gihawi, A., Cardenas, R., Hurst, R., & Brewer, D. S. (2023). Quality Control in Metagenomics Data. In S. Mitra (Ed.), Metagenomic Data Analysis (pp. 21–54). Springer US. https://doi.org/10.1007/978-1-0716-3072-3_2

[5] Chen, S. (2025). fastp 1.0: An ultra‐fast all‐round tool for FASTQ data quality control and preprocessing. iMeta, 4(5), e70078-n/a. https://doi.org/10.1002/imt2.70078

[6] Wood, D. E., Lu, J., & Langmead, B. (2019). Improved metagenomic analysis with Kraken 2. Genome Biology, 20(1), Article 257. https://doi.org/10.1186/s13059-019-1891-0

[7] Wright, R. J., Comeau, A. M., & Langille, M. G. I. (2023). From defaults to databases: parameter and database choice dramatically impact the performance of metagenomic taxonomic classification tools. Microbial Genomics, 9(3), Article 000949. https://doi.org/10.1099/mgen.0.000949

[8] Lu, J., Breitwieser, F. P., Thielen, P., & Salzberg, S. L. (2017). Bracken: estimating species abundance in metagenomics data. PeerJ. Computer Science, 3, Article e104. https://doi.org/10.7717/peerj-cs.104

[9] McMURDIE, P. J., & Holmes, S. (2012). Phyloseq: a bioconductor package for handling and analysis of high-throughput phylogenetic sequence data. In Biocomputing 2012 (pp. 235-246).

[10] Lin, H., & Peddada, S. D. (2024). Multigroup analysis of compositions of microbiomes with covariate adjustments and repeated measures. Nature Methods, 21(1), 83–91. https://doi.org/10.1038/s41592-023-02092-7

[11] Ncbi. (n.d.). NCBI/SRA-Tools: Sra tools. GitHub. https://github.com/ncbi/sra-tools

[12] Ncbi. (n.d.). 08. Prefetch and fasterq dump. GitHub. https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump 

[13] S-Andrews. (n.d.). S-andrews/FASTQC: A quality control analysis tool for high throughput sequencing data. GitHub. https://github.com/s-andrews/fastqc 

[14] OpenGene. (n.d.). OpenGene/fastp: An ultra-fast all-in-one FASTQ preprocessor (QC/adapters/trimming/filtering/splitting/merging...). GitHub. https://github.com/opengene/fastp 

[15] DerrickWood. (n.d.). Derrickwood/Kraken2: The second version of the Kraken Taxonomic Sequence Classification System. GitHub. https://github.com/DerrickWood/kraken2 

[16] Index zone. Index zone by BenLangmead. (n.d.). https://benlangmead.github.io/aws-indexes/k2 

[17] jenniferlu717. (n.d.). Jenniferlu717/Bracken: Bracken (bayesian reestimation of abundance with Kraken) is a highly accurate statistical method that computes the abundance of species in DNA sequences from a metagenomics sample. GitHub. https://github.com/jenniferlu717/Bracken 

[18] joey711. (n.d.). Joey711/phyloseq: Phyloseq is a set of classes, wrappers, and tools (in R) to make it easier to import, store, and analyze phylogenetic sequencing data; and to reproducibly share that data and analysis with others. see the phyloseq front page:. GitHub. https://github.com/joey711/phyloseq 

[19] Cran. (n.d.). Cran/vegan: :exclamation: This is a read-only mirror of the cran R package repository. vegan - community ecology package. homepage: Https://vegandevs.github.io/vegan/, https://github.com/vegandevs/vegan report bugs for this package: Https://github.com/vegandevs/vegan/issues. GitHub. https://github.com/cran/vegan/tree/master 

[20] Tidyverse. (n.d.). Tidyverse/GGPLOT2: An implementation of the grammar of graphics in R. GitHub. https://github.com/tidyverse/ggplot2 

[21] FrederickHuangLin. (n.d.). Frederickhuanglin/ANCOMBC: Differential Abundance (DA) and correlation analyses for microbial absolute abundance data. GitHub. https://github.com/FrederickHuangLin/ANCOMBC 
