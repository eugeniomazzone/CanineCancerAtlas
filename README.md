# The Canine Cancer Genome Atlas

## DISCLAIMER
All scripts presented here were tested to work properly for PRJNA752630 BioProject. Before applying the code to other datasets, carefully go through the scrips to make sure everythig is correctly set.
This applies specifically to the **_RunPreProcess.sh_** script where SRA metadata are accessed (SraRunTable.txt). 
Moreover, to perform statistical analyses on two (or more) histotypes/BioPorjects simultaneusly, please append them in the R scripts headers (reads as the following).
```
biop <- c('PRJNA752630')
tt <- c('B-cell Lymphoma')
```

## General Informations

The Canine Cancer Genome Atlas (CCGA) is a project that aims to build and mantain a comprehensive repository of genetic aberrations found in canine neoplasms.
This repository aims to give an overview of the bioinformatic tools and detail the core statistical analysis performed on the data.
To have an in depth explanation of the method, please refer to [The Genetic Landscape of Canine Tumors: Insights from the Canine Cancer Genome Atlas (CCGA)](https://doi.org/10.21203/rs.3.rs-5025541/v1).

## Data selection and inclusion criteria

The following pipeline is intended to be used for **canine cancer data** produced in **Next Generation Sequencing** (NGS) experiment.
Moreover, all the data selected for the study where included if and only if they were complaiant with all the following specifications:
- Sequencing data from healty tissue from the same animal were available (mathced normal).
- Tumor mean coverage was at least 15x.
- Matched normal covearge was at least 10x.

## Data Preprocessing

Data preprocessing follow the Genome Analysis Toolkit (GATK) guidelines. Briefly:
- Data are aquired from the Sequence Read Archive (SRA).
- Converted to FASTQ format (sra-toolkit).
- Aligned to the reference genome using Burrow-Wheeler Aligner (BWA).
- Converted to BAM format.
- Sorted by Coordinate (with Samtools).
- Read groups information added.
- The duplicates were marked (with Picard).
- Base Quality Score Recalibration (BQRS) was performed.

## Variant Calling

Our pipeline serves to call three variant categories: SNV and InDels, Copy-Number Aberrations (CNA) and Structural Variants (SV).

1. Single-Nucleotide Variants (SNV) and small insertions and deletions (InDel) were called simultaneusly with three different tools, Mutect,Strelka and Varscan. Then the results were merged using a classic majority voting procedure where only variables found by at least two-out-of-three caller were kept for downstream analysis.

2. ASCAT, FACETS and SEQUENZA were used to access copy number information for each samples. They were installed and used as part of the EaCoN R packege.

3. Finally, structural vaiants were retrived with Delly program (dellytools).

## Annotations and Statistical Analysis

Annotations were performed using ANNOVAR (SNV and InDel), GTF files plus R scripts (CNA) and sansa (SV).
Statistical analysis were performed to access many aspects of the data.

## Additional resources to run the analysis

To run the many scripts above you'll need to download locally the following files:
- CanFam3.1 reference genome.
- CanFam3.1 cytoband information.
- CanFam3.1 Gene Name GTF.
- Germline SNP data from DogSD website.

## Website

 For an in-depth (and user friendly) exploration the results, visit [The Canine Cancer Genome Atlas website](https://caninecancergenomeatlas.org/).

## Citing

If you reuse code from this repository, please cite [The Genetic Landscape of Canine Tumors: Insights from the Canine Cancer Genome Atlas (CCGA)](https://doi.org/10.21203/rs.3.rs-5025541/v1).

