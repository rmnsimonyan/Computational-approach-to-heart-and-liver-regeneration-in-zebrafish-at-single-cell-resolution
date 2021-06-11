library(Matrix)
library(Seurat)
library(DESeq2)

#Read the data
dat <- Read10X(data.dir = "ZF/HEP/")
colnames(dat) = paste0("HEP", colnames(dat) )
HEP <- CreateSeuratObject(counts = dat, min.cells = 5, min.features = 200, project = "HEP")

dat2 <- Read10X(data.dir = "ZF/MTZablhep/")
colnames(dat2) = paste0("MTZabl", colnames(dat2) )
MTZabl <- CreateSeuratObject(counts = dat2, min.cells = 5, min.features = 200, project = "MTZabl")


# Check the percentage of mt genes
HEP[["percent.mt"]] <- PercentageFeatureSet(HEP, pattern = "^mt-")
VlnPlot(HEP, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

MTZabl[["percent.mt"]] <- PercentageFeatureSet(MTZabl, pattern = "^mt-")
VlnPlot(MTZabl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


#Clustering

# QC
qHEP <- subset(HEP, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
qMTZabl <- subset(MTZabl, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

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

qHEP <- JackStraw(qHEP, num.replicate = 100)
qHEP <- ScoreJackStraw(qHEP, dims = 1:20)
JackStrawPlot(qHEP, dims = 1:15)

ElbowPlot(qHEP)


#Cluster

#Control
qHEP <- FindNeighbors(qHEP, dims = 1:10)
qHEP <- FindClusters(qHEP, resolution = 0.5)
qHEP <- RunUMAP(qHEP, dims = 1:10)
DimPlot(qHEP, reduction = "umap")

cluster1.markers <- FindMarkers(qHEP, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 20)

cluster2.markers <- FindMarkers(qHEP, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 20)

cluster4.markers <- FindMarkers(qHEP, ident.1 = 4, min.pct = 0.25)
head(cluster4.markers, n = 20)

VlnPlot(qHEP, features = c("fabp10a", "selenop", "hnf4", "gc", "cp"))

VlnPlot(qHEP, features = c("anxa4","epcam", "sox9b"))


#Injured
qMTZabl <- FindNeighbors(qMTZabl, dims = 1:10)
qMTZabl <- FindClusters(qMTZabl, resolution = 0.5)
qMTZabl <- RunUMAP(qMTZabl, dims = 1:10)
DimPlot(qMTZabl, reduction = "umap")

cluster2.markers <- FindMarkers(qMTZabl, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 10)

cluster4.markers <- FindMarkers(qMTZabl, ident.1 = 4, min.pct = 0.25)
head(cluster4.markers, n = 10)

VlnPlot(qMTZabl, features = c("fabp10a", "selenop", "hnf4a", "gc", "cp", "anxa4","epcam", "sox9b"))

VlnPlot(qMTZabl, features = c("anxa4","epcam", "sox9b"))

VlnPlot(qHEP, features = c("cp1","hmgcra", "fads2"))

#Farnsworth 2019 cluster 55
VlnPlot(qHEP, features = c("cbln14",	"si:ch211-113d11.5",	"zgc:171534",	"or109-11",	"serpina7",	"cbln5",	"cbln9",	"ccl39.2",	"zgc:172253",	"cpb2",	"ahsg1",	"si:ch211-262h13.5",	"si:ch73-15n24.1",	"si:dkey-96f10.1",	"itih2",	"si:dkeyp-73d8.6"))
VlnPlot(qMTZabl, features = c("cbln14",	"si:ch211-113d11.5",	"zgc:171534",	"or109-11",	"serpina7",	"cbln5",	"cbln9",	"ccl39.2",	"zgc:172253",	"cpb2",	"ahsg1",	"si:ch211-262h13.5",	"si:ch73-15n24.1",	"si:dkey-96f10.1",	"itih2",	"si:dkeyp-73d8.6"))


#Farnsworth 2019 cluster 121
VlnPlot(qHEP, features = c("si:ch211-186e20.2",	"cbln17",	"cyp2x10.2",	"tdrd5",	"apobb.2",	"AL953907.1",	"crp2",	"igfals",	"zgc:103559",	"cfhl3",	"uroc1",	"gys2",	"f9b",	"acsbg1",	"CR626907.1",	"fgf19"))
VlnPlot(qMTZabl, features = c("si:ch211-186e20.2",	"cbln17",	"cyp2x10.2",	"tdrd5",	"apobb.2",	"AL953907.1",	"crp2",	"igfals",	"zgc:103559",	"cfhl3",	"uroc1",	"gys2",	"f9b",	"acsbg1",	"CR626907.1",	"fgf19"))

#Farnsworth 2019 cluster 217
VlnPlot(qHEP, features = c("si:dkey-79f11.10",	"BX927210.1",	"CABZ01000634.1",	"si:ch211-113d11.8",	"si:dkey-7f3.14",	"FO834888.1",	"CU856343.1",	"gltpd2",	"rln3a",	"acod1",	"soat2",	"apoc4",	"gpr84",	"hp",	"si:ch211-271e10.2",	"c3b.1"))
VlnPlot(qMTZabl, features = c("si:dkey-79f11.10",	"BX927210.1",	"CABZ01000634.1",	"si:ch211-113d11.8",	"si:dkey-7f3.14",	"FO834888.1",	"CU856343.1",	"gltpd2",	"rln3a",	"acod1",	"soat2",	"apoc4",	"gpr84",	"hp",	"si:ch211-271e10.2",	"c3b.1"))

VlnPlot(qMTZabl,features = c("cyp3c1", "mmp7"))

VlnPlot(qMTZabl, features = c("glula",
                           "cyp2x8",
                           "ass1",
                           "asl",
                           "axin2",
                           "zgc:173994",
                           "cyp1d1",
                           "psmd4a",
                           "arg1",
                           "pck1",
                           "sdhda"))

VlnPlot(qHEP, features = c("atf5a",
                           "hamp",
                           "fgl1",
                           "fabp1a",
                           "serpina1",
                           "rbp4",
                           "ass1",
                           "agxta",
                           "rnasel3",
                           "hpxa",
                           "hpdb",
                           "fga",
                           "cyp1d1",
                           "cyp4f3",
                           "apoba",
                           "mlxipl",
                           "cp",
                           "cyp2aa8",
                           "crp",
                           "apoea",
                           "serpinf2a",
                           "c1r",
                           "asgrl1",
                           "itih3a.1"))




