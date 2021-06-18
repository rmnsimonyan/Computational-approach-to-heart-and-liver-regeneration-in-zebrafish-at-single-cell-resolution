# Computational approach to heart and liver regeneration in zebrafish at single cell resolution

This project aimed to run several clustering algorithms on healthy and ablated liver scRNA-seq data from 10XGenomics. "data" folder contains compressed sample scRNA-seq data from 10XGenomics that could be used to run the pipeline on. In the same folder there is a sample expression matrix which can be obtained from 10XGenomics format. The Seurat pipeline runs on the 10XGenomics format data (that needs to be decompressed in the same folder before running). Meanwhile, scPopCorn and SOM pipelines run on the expression matrix txt files (data1.txt data2.txt).

How to run:

## Seurate clustering script

1. Decompressing the data files

```
cd data
unzip T.zip
unzip NK.zip
```

2. Run Seurate pipeline that does standard and integrated clustering if the datasets and writes the results in the files (graphs, tables). The enrichment (ORA) results of the provided sample dataset are not significant, that's why it is not going to provide any enrichment results if run on the provided dataset. However, on the actual dataset the same ORA pipeline works.

```
Rscript Seurat_analysis.R
```
## scPopcorn clustering

```
python scPopcorn.py
```

## SOM clustering

```
 Rscipt SOM.R
```
