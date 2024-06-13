# MSc_proteins-crc

## overview
There is growing evidence, from Mendelian randomisation analyses, that some proteins are "causal" for colorectal cancer (CRC) development. A limitation of Mendelian randomisation analyses is that the data used (plasma proteins) do not necessarily reflect a localised exposure, for instance protein abundance measured in serum may not reflect protein abundance in the colon. Tissue-specific MR uses tissue specific (e.g., colon tissue) gene expression data to identify instruments to be used in MR analyses. This project aims to establish whether plasma proteins causal for CRC are associated with colorectal cancer in colon tissue.

## data
We obtained data for (1) plasma proteins, (2) gene expression, and (3) CRC from:

1. plasma proteins are from [Ferkingstad et al., (2021)](https://www.nature.com/articles/s41588-021-00978-w) and were downloaded from [deCODE](https://download.decode.is/form/folder/proteomics)
    * we used the protein [GREM1](https://www.ncbi.nlm.nih.gov/gene/26585) as an exemplar

2. gene expression data are from [GTEx V8](https://www.gtexportal.org/home/downloads/adult-gtex/qtl)
    * [README](https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl/README_eQTL_v8.txt)
    * we used the protein [GREM1](https://gtexportal.org/home/gene/GREM1) as an exemplar
    * we used data from colon tissue (sigmoid and transverse)

3. CRC data are from [Fernandez-Rozadilla et al., 2023](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10094749/) and were downloaded from the [GWAS catalog](https://www.ebi.ac.uk/gwas/studies/GCST90129505) (accession: GCST90129505)

## analysis
1. MR analysis of plasma proteins and CRC
2. MR analysis of CRC and plasma proteins
2. colocalisation of plasma proteins and CRC
3. colocalisation of plasma proteins and gene expression
4. colocalisation of gene expression and plasma proteins
5. MR analysis of gene expression and CRC
6. Gene-Set enrichment analysis of associated proteins
7. Sensitivity analysis 

## code

