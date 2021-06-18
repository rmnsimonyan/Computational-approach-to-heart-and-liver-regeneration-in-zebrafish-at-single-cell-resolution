library(Matrix)
library(Seurat)
library(DESeq2)

dir.create("Seurat_standard")

#Read the data
dat <- Read10X(data.dir = "T/")
colnames(dat) = paste0("T", colnames(dat) )
HEP <- CreateSeuratObject(counts = dat, min.cells = 5, min.features = 200, project = "T")

dat2 <- Read10X(data.dir = "NK/")
colnames(dat2) = paste0("NK", colnames(dat2) )
MTZabl <- CreateSeuratObject(counts = dat2, min.cells = 5, min.features = 200, project = "NK")


# Check the percentage of mt genes
HEP[["percent.mt"]] <- PercentageFeatureSet(HEP, pattern = "^MT-")

jpeg("Seurat_standard/Seurat_dat1_stats.jpeg")
VlnPlot(HEP, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()


MTZabl[["percent.mt"]] <- PercentageFeatureSet(MTZabl, pattern = "^MT-")

jpeg("Seurat_standard/Seurat_dat2_stats.jpeg")
VlnPlot(MTZabl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

#Clustering

# QC
qHEP <- subset(HEP, subset = nFeature_RNA > 200)
qMTZabl <- subset(MTZabl, subset = nFeature_RNA > 200) 

# For filtering hepatocytes
#qHEP = subset(qHEP, subset = fabp10a > 1)
#qMTZabl = subset(qMTZabl, subset = fabp10a > 1)

#Normalize abd scale
qHEP <- NormalizeData(qHEP)
qMTZabl <- NormalizeData(qMTZabl)

qHEP <- ScaleData(qHEP, features = rownames(qHEP))
qMTZabl <- ScaleData(qMTZabl, features = rownames(qMTZabl))

#Find highly variable genes
qHEP <- FindVariableFeatures(qHEP)
qHEP <- RunPCA(qHEP, features = VariableFeatures(object = qHEP))

qMTZabl <- FindVariableFeatures(qMTZabl)
qMTZabl <- RunPCA(qMTZabl, features = VariableFeatures(object = qMTZabl))


VizDimLoadings(qHEP, dims = 1:2, reduction = "pca")
DimPlot(qHEP, reduction = "pca")

VizDimLoadings(qMTZabl, dims = 1:2, reduction = "pca")
DimPlot(qMTZabl, reduction = "pca")

# Checking dimentionality of datasets
#qHEP <- JackStraw(qHEP, num.replicate = 100)
#qHEP <- ScoreJackStraw(qHEP, dims = 1:20)
#JackStrawPlot(qHEP, dims = 1:15)

#ElbowPlot(qHEP)


#Cluster

#Control
qHEP <- FindNeighbors(qHEP, dims = 1:10)
qHEP <- FindClusters(qHEP, resolution = 0.5)
qHEP <- RunUMAP(qHEP, dims = 1:10)

jpeg("Seurat_standard/Seurat_dat1_clustering.jpeg")
DimPlot(qHEP, reduction = "umap")
dev.off()

cluster0.markers <- FindMarkers(qHEP, ident.1 = 0, min.pct = 0.25)
head(cluster0.markers, n = 20)
write.csv(cluster0.markers, file="Seurat_dat1_cluster0.csv")

cluster1.markers <- FindMarkers(qHEP, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 20)
write.csv(cluster1.markers, file="Seurat_dat1_cluster1.csv")

cluster2.markers <- FindMarkers(qHEP, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 20)
write.csv(cluster2.markers, file="Seurat_dat1_cluster2.csv")

cluster3.markers <- FindMarkers(qHEP, ident.1 = 3, min.pct = 0.25)
head(cluster3.markers, n = 20)
write.csv(cluster3.markers, file="Seurat_dat1_cluster3.csv")

cluster4.markers <- FindMarkers(qHEP, ident.1 = 4, min.pct = 0.25)
head(cluster4.markers, n = 20)
write.csv(cluster4.markers, file="Seurat_dat1_cluster4.csv")

cluster5.markers <- FindMarkers(qHEP, ident.1 = 5, min.pct = 0.25)
head(cluster5.markers, n = 20)
write.csv(cluster5.markers, file="Seurat_dat1_cluster5.csv")

dat1clusters = list(cluster0.markers, cluster1.markers, cluster2.markers, cluster3.markers,
                    cluster4.markers, cluster5.markers)

# Check some marker genes' expression  patterns across clusters

jpeg("Seurat_standard/Seurat_dat1_markergenes_in_clusters.jpeg")
VlnPlot(qHEP, features = c("HLA-B","AIF1", "CCL5"))
dev.off()

# marker genes for zebrafish heptocytes and BECs 
#VlnPlot(qHEP, features = c("fabp10a", "selenop", "hnf4", "gc", "cp"))
#VlnPlot(qHEP, features = c("anxa4","epcam", "sox9b"))


#Injured
qMTZabl <- FindNeighbors(qMTZabl, dims = 1:10)
qMTZabl <- FindClusters(qMTZabl, resolution = 0.5)
qMTZabl <- RunUMAP(qMTZabl, dims = 1:10)

jpeg("Seurat_standard/Seurat_dat2_clustering.jpeg")
DimPlot(qMTZabl, reduction = "umap")
dev.off()


cluster0.markers <- FindMarkers(qMTZabl, ident.1 = 0, min.pct = 0.25)
head(cluster0.markers, n = 20)
write.csv(cluster0.markers, file="Seurat_standard/Seurat_dat2_cluster0.csv")

cluster1.markers <- FindMarkers(qMTZabl, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 20)
write.csv(cluster1.markers, file="Seurat_standard/Seurat_dat2_cluster1.csv")

cluster2.markers <- FindMarkers(qMTZabl, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 20)
write.csv(cluster2.markers, file="Seurat_standard/Seurat_dat2_cluster2.csv")

cluster3.markers <- FindMarkers(qMTZabl, ident.1 = 3, min.pct = 0.25)
head(cluster3.markers, n = 20)
write.csv(cluster3.markers, file="Seurat_standard/Seurat_dat2_cluster3.csv")

cluster4.markers <- FindMarkers(qMTZabl, ident.1 = 4, min.pct = 0.25)
head(cluster4.markers, n = 20)
write.csv(cluster4.markers, file="Seurat_standard/Seurat_dat2_cluster4.csv")

cluster5.markers <- FindMarkers(qMTZabl, ident.1 = 5, min.pct = 0.25)
head(cluster5.markers, n = 20)
write.csv(cluster5.markers, file="Seurat_standard/Seurat_dat2_cluster5.csv")

cluster6.markers <- FindMarkers(qMTZabl, ident.1 = 6, min.pct = 0.25)
head(cluster6.markers, n = 20)
write.csv(cluster6.markers, file="Seurat_standard/Seurat_dat2_cluster6.csv")

dat2clusters = list(cluster0.markers, cluster1.markers, cluster2.markers, cluster3.markers,
                    cluster4.markers, cluster5.markers, cluster6.markers)
# Check some marker genes' expression  patterns across clusters

jpeg("Seurat_standard/Seurat_dat2_markergenes_in_clusters.jpeg")
VlnPlot(qMTZabl, features = c("HLA-B","AIF1", "CCL5"))
dev.off()

#VlnPlot(qMTZabl, features = c("fabp10a", "selenop", "hnf4a", "gc", "cp", "anxa4","epcam", "sox9b"))
#VlnPlot(qMTZabl, features = c("anxa4","epcam", "sox9b"))
#VlnPlot(qHEP, features = c("cp1","hmgcra", "fads2"))

#Farnsworth 2019 cluster 55
#VlnPlot(qHEP, features = c("cbln14",	"si:ch211-113d11.5",	"zgc:171534",	"or109-11",	"serpina7",	"cbln5",	"cbln9",	"ccl39.2",	"zgc:172253",	"cpb2",	"ahsg1",	"si:ch211-262h13.5",	"si:ch73-15n24.1",	"si:dkey-96f10.1",	"itih2",	"si:dkeyp-73d8.6"))
#VlnPlot(qMTZabl, features = c("cbln14",	"si:ch211-113d11.5",	"zgc:171534",	"or109-11",	"serpina7",	"cbln5",	"cbln9",	"ccl39.2",	"zgc:172253",	"cpb2",	"ahsg1",	"si:ch211-262h13.5",	"si:ch73-15n24.1",	"si:dkey-96f10.1",	"itih2",	"si:dkeyp-73d8.6"))


#Farnsworth 2019 cluster 121
#VlnPlot(qHEP, features = c("si:ch211-186e20.2",	"cbln17",	"cyp2x10.2",	"tdrd5",	"apobb.2",	"AL953907.1",	"crp2",	"igfals",	"zgc:103559",	"cfhl3",	"uroc1",	"gys2",	"f9b",	"acsbg1",	"CR626907.1",	"fgf19"))
#VlnPlot(qMTZabl, features = c("si:ch211-186e20.2",	"cbln17",	"cyp2x10.2",	"tdrd5",	"apobb.2",	"AL953907.1",	"crp2",	"igfals",	"zgc:103559",	"cfhl3",	"uroc1",	"gys2",	"f9b",	"acsbg1",	"CR626907.1",	"fgf19"))

#Farnsworth 2019 cluster 217
#VlnPlot(qHEP, features = c("si:dkey-79f11.10",	"BX927210.1",	"CABZ01000634.1",	"si:ch211-113d11.8",	"si:dkey-7f3.14",	"FO834888.1",	"CU856343.1",	"gltpd2",	"rln3a",	"acod1",	"soat2",	"apoc4",	"gpr84",	"hp",	"si:ch211-271e10.2",	"c3b.1"))
#VlnPlot(qMTZabl, features = c("si:dkey-79f11.10",	"BX927210.1",	"CABZ01000634.1",	"si:ch211-113d11.8",	"si:dkey-7f3.14",	"FO834888.1",	"CU856343.1",	"gltpd2",	"rln3a",	"acod1",	"soat2",	"apoc4",	"gpr84",	"hp",	"si:ch211-271e10.2",	"c3b.1"))

#VlnPlot(qMTZabl,features = c("cyp3c1", "mmp7"))

#VlnPlot(qMTZabl, features = c("glula",
                              # "cyp2x8",
                              # "ass1",
                              # "asl",
                              # "axin2",
                              # "zgc:173994",
                              # "cyp1d1",
                              # "psmd4a",
                              # "arg1",
                              # "pck1",
                              # "sdhda"))
# 
# VlnPlot(qHEP, features = c("atf5a",
#                            "hamp",
#                            "fgl1",
#                            "fabp1a",
#                            "serpina1",
#                            "rbp4",
#                            "ass1",
#                            "agxta",
#                            "rnasel3",
#                            "hpxa",
#                            "hpdb",
#                            "fga",
#                            "cyp1d1",
#                            "cyp4f3",
#                            "apoba",
#                            "mlxipl",
#                            "cp",
#                            "cyp2aa8",
#                            "crp",
#                            "apoea",
#                            "serpinf2a",
#                            "c1r",
#                            "asgrl1",
#                            "itih3a.1"))
# 
# #mitotic cell cycle genes from GOterm from GSEA analysis results
# VlnPlot(qMTZabl, features = c("btg1",
#                               "tubb4b",
#                               "vcp",
#                               "ran",
#                               "tuba8l2",
#                               "si:dkey-79d12.5",
#                               "chmp4ba",
#                               "btg2",
#                               "ccni"))
# 
# 

# Over-representation Analysis for cluster 6
# organism = "org.Dm.eg.db"
organism = "org.Hs.eg.db"

#BiocManager::install(organism, character.only = TRUE)

library(organism, character.only = TRUE)
library(clusterProfiler)

for(i in dat1clusters){
  go_enrich <- enrichGO(gene = rownames(i)[i$p_val_adj<0.05],
                        OrgDb = organism, 
                        keyType = "GENENAME",
                        readable = T,
                        ont = "ALL",
                        pvalueCutoff = 0.05, 
                        qvalueCutoff = 0.10)
  print(go_enrich)
  write.table(go_enrich[,3], paste("Seurat_standard/Seurat_dat1_Enrichment_clust_",i,".txt",sep=""))
  
}


for(i in dat2clusters){
go_enrich <- enrichGO(gene = rownames(i)[i$p_val_adj<0.05],
                       OrgDb = organism, 
                       keyType = "GENENAME",
                       readable = T,
                       ont = "ALL",
                       pvalueCutoff = 0.05, 
                       qvalueCutoff = 0.10)
print(go_enrich)
write.table(go_enrich[,3], paste("Seurat_standard/Seurat_dat2_Enrichment_clust_",i,".txt",sep=""))

}


# Integrated clustering
# install.packages("https://seurat.nygenome.org/src/contrib/ifnb.SeuratData_3.0.0.tar.gz", repos = NULL, type = "source") 
library(SeuratData)
library(patchwork)
library(ifnb.SeuratData)

dir.create("Seurat_Integrated")

# Set up control object
dat <- Read10X(data.dir = "T/")
colnames(dat) = paste0("T", colnames(dat) )
HEP <- CreateSeuratObject(counts = dat, min.cells = 5, min.features = 200, project = "T")
HEPq <- subset(HEP, subset = nFeature_RNA > 200) #& fabp10a > 1)
HEPq$stim <- "HEP"
HEPq <- NormalizeData(HEPq, verbose = FALSE)
HEPq <- FindVariableFeatures(HEPq, selection.method = "vst", nfeatures = 2000)

# Set up stimulated object
dat2 <- Read10X(data.dir = "NK/")
colnames(dat2) = paste0("NK", colnames(dat2) )
MTZabl <- CreateSeuratObject(counts = dat2, min.cells = 5, min.features = 200, project = "NK")
MTZablq <- subset(MTZabl, subset = nFeature_RNA > 200) #& fabp10a > 1)
MTZablq$stim <- "MTZablq"
MTZablq <- NormalizeData(MTZablq, verbose = FALSE)
MTZablq <- FindVariableFeatures(MTZablq, selection.method = "vst", nfeatures = 2000)

immune.anchors <- FindIntegrationAnchors(object.list = list(HEPq, MTZablq), dims = 1:20)

immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)

# specify that we will perform downstream analysis on the corrected data note that the original
# unmodified data still resides in the 'RNA' assay
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

# Visualization

jpeg("Seurat_Integrated/Seurat_Integrated_clustering.jpeg")
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
dev.off()

jpeg("Seurat_Integrated/Seurat_Integrated_clustering_seprately.jpeg")
DimPlot(immune.combined, reduction = "umap", split.by = "stim")
dev.off()


organism = "org.Hs.eg.db"
#BiocManager::install(organism, character.only = TRUE)

library(organism, character.only = TRUE)
library(clusterProfiler)

markers = list()
enriched = list()

for(i in 0:(length(unique(Idents(immune.combined)))-1)){
  markers[[i]] = FindMarkers(immune.combined, ident.1 = i, min.pct = 0.25)
  head(markers[[i]], n = 20)
  write.csv(markers[[i]], file=paste("Seurat_Integrated/Seurat_Integrated_cluster_",i,".csv", sep = ""))
  
  go_enrich <- enrichGO(gene = rownames(markers[[i]])[markers[[i]]$p_val_adj<0.05],
                        #OrgDb = organism, 
                        keyType = "GENENAME",
                        readable = T,
                        ont = "ALL",
                        pvalueCutoff = 0.05, 
                        qvalueCutoff = 0.10)
  print(go_enrich[,2])
  write.table(go_enrich[,3], paste("Seurat_Integrated/Seurat_integrated_Enrichment_clust_",i,".txt",sep=""))
  }
