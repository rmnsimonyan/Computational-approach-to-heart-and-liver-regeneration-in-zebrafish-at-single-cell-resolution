# Self organizing map clustering with oposSOM package

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#
#BiocManager::install("oposSOM")

rm(list = ls())
library(oposSOM)
library(stringr)

HEPfin = read.delim("data/data1.txt", header = TRUE, row.names = 1)

MTZablfin = read.delim("data/data2.txt", header = TRUE, row.names = 1)


mat = merge(HEPfin,MTZablfin,by="row.names",all=TRUE)
rownames(mat) = mat[,1]
mat = mat[,-1]
mat[is.na(mat)] = 0

#mat = read.csv("exptable.csv", header = TRUE, row.names = 1)

genenumber = c()
for(i in colnames(mat)){
  genenumber[i] = sum(mat[,i]>0)
}


mat = mat[,names(genenumber[which(!genenumber<200)])]

env <- opossom.new(list(dataset.name = "SD.som.maf.Diagnosis_new.onset",
                        
                        dim.1stLvlSom = "auto",
                        dim.2ndLvlSom = 10,
                        
                        training.extension = 1,
                        rotate.SOM.portraits = 0,
                        flip.SOM.portraits = F,
                        
                        #database.biomart = "drerio_gene_ensembl",# "ENSEMBL_MART_ENSEMBL",
                        #database.host = "" ,#"jul2015.archive.ensembl.org",
                        database.dataset = "drerio_gene_ensembl",  #"auto",
                        #database.id.type = "",
                        
                        activated.modules = list( "reporting" = T,
                                                  "primary.analysis" = T, 
                                                  "sample.similarity.analysis" = T,
                                                  "geneset.analysis" = F, 
                                                  "geneset.analysis.exact" = F,
                                                  "group.analysis" = T,
                                                  "difference.analysis" = T,
                                                  "psf.analysis" = F ),
                        
                        standard.spot.modules = "group.overexpression",
                        
                        spot.coresize.modules = 7,
                        spot.threshold.modules = 0.95,
                        spot.coresize.groupmap = 10,
                        spot.threshold.groupmap = 0.90,
                        
                        feature.centralization = T,
                        sample.quantile.normalization = T,
                        
                        pairwise.comparison.list = list() ) )

env$indata = mat

labels = c(rep("HEP",sum(str_detect(colnames(mat), "HEP")) ), rep("MTZabl",sum(str_detect(colnames(mat), "MTZabl")) ))
names(labels) = colnames(mat)
env$group.labels <- labels

opossom.run(env)

#Enrich one of the group overexpression spots as an example

organism = "org.Hs.eg.db"
#BiocManager::install(organism, character.only = TRUE)

library(organism, character.only = TRUE)
library(clusterProfiler)

go_enrichE <- enrichGO(gene = env$spot.list.group.overexpression$spots$A$genes,
                       OrgDb = organism, 
                       keyType = "GENENAME",
                       readable = T,
                       ont = "ALL",
                       pvalueCutoff = 0.05, 
                       qvalueCutoff = 0.10)
print(go_enrichE[,3])
