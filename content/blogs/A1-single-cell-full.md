---
title: 一文读懂单细胞测序分析流程
date: 16th Apr 2019
description: 单细胞测序中文教程.
image: /blogs-img/singlecell.webp
alt: 一文读懂单细胞测序分析流程
ogImage: /blog-img/singlecell.webp
tags: ['单细胞测序', '生物信息', 'r语言']
published: true
---

- [x] 本教程已更新，更新时间：
- [x] 2019/12/27

### 摘要

一文介绍单细胞测序生物信息分析完整流程，这可能是最新也是最全的流程

### 基础流程（cellranger）

![](https://i.loli.net/2019/12/27/wMSmukTrLyaBeql.png)

### cellranger 数据拆分

`cellranger mkfastq` 可用于将单细胞测序获得的 BCL 文件拆分为可以识别的 fastq 测序数据

```bash
cellranger makefastq   --run=[ ]   --samplesheet=[sample.csv] --jobmode=local --localcores=20 --localmem=80
```

> -–run ：是下机数据 BCL 所在的路径；
> -–samplesheet ：样品信息列表--共三列（lane id ,sample name ,index name)
> 注意要控制好核心数和内存数

运行产出结果存在于 out 目录中

### cellranger 数据统计

`cellranger count`是 cellranger 最主要也是最重要的功能：完成细胞和基因的定量，也就是产生了我们用来做各种分析的基因表达矩阵。

```bash
cellranger count \
-–id=sample345 \
-–transcriptome=/opt/refdata-cellranger-GRCh38-1.2.0/GRCh38 \
-–fastqs=/home/jdoe/runs/HAWT7ADXX/outs/fastq_path \
-–indices=SI-3A-A1 \
–-cells=1000
```

> id ：产生的结果都在这个文件中，可以取几号样品（如 sample345）；

> fastqs ：由 cellranger mkfastq 产生的 fastqs 文件夹所在的路径；fastqs ：由 cellranger mkfastq 产生的 fastqs 文件夹所在的路径；

> indices：sample index：SI-3A-A1；

> transcriptome：参考转录组文件路径；

> cells：预期回复的细胞数；

### 下游分析

cellranger count 计算的结果只能作为初步观测的结果，如果需要进一步分析聚类细胞，还需要进行下游分析，这里使用官方推荐 R 包（Seurat 3.1）

流程参考官方（[外周血分析标准流程](https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html)）

### 软件安装

```bash
install.packages('Seurat')
library(Seurat)
```

### 生成 Seruat 对象

```r
library(dplyr)
library(Seurat)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "../data/pbmc3k/filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc
```

```bash
## An object of class Seurat
## 13714 features across 2700 samples within 1 assay
## Active assay: RNA (13714 features)
```

> 这里读取的是单细胞 count 结果中的矩阵目录；
> 在对象生成的过程中，做了初步的过滤；
> 留下所有在>=3 个细胞中表达的基因 min.cells = 3；
> 为了除去一些质量差的细胞,留下所有检测到>=200 个基因的细胞 min.genes = 200。

### 标准预处理流程

```r
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")
```

这一步 mit-开头的为线粒体基因，这里将其进行标记并统计其分布频率

```r
# Visualize QC metrics as a violin plot
VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

对 pbmc 对象做小提琴图，分别为基因数，细胞数和线粒体占比

![](https://i.loli.net/2019/12/27/KgpU8PVEWylzZAx.png)

```r
pbmc <- subset(x = pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
```

接下来，根据图片中基因数和线粒体数，分别设置过滤参数，这里基因数 200-2500，线粒体百分比为小于 5%

### 数据标准化

```r
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- NormalizeData(object = pbmc)
```

### 鉴定高度变化基因

```r
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(x = VariableFeatures(object = pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(object = pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
```

![](https://i.loli.net/2019/12/27/U31g6tmwEfaYX78.png)

### 数据归一化

```r
all.genes <- rownames(x = pbmc)
pbmc <- ScaleData(object = pbmc, features = all.genes)
```

这里设置对所有的基因都做了`scale`,但是需要知道的是，其实后续的分析都是基于高变基因的，因此，使用默认参数就可以了，而且提升效率。

```r
pbmc <- ScaleData(object = pbmc)
```

### 线形降维

```r
pbmc <- RunPCA(object = pbmc, features = VariableFeatures(object = pbmc))
```

这里有多种方法展示 pca 结果，本文采用最简单的方法

```r
DimPlot(object = pbmc, reduction = "pca")
```

![](https://i.loli.net/2019/12/27/UyOzTfhNFBMt4rg.png)

### 鉴定数据集的可用维度

```r
pbmc <- JackStraw(object = pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(object = pbmc, dims = 1:20)
JackStrawPlot(object = pbmc, dims = 1:15)
```

![](https://i.loli.net/2019/12/27/3xBbj7RHtPsyniK.png)

虚线以上的为可用维度，你也可以调整 dims 参数，画出所有 pca 查看

另外一种鉴定手段是绘制所有 PC 的分布点图

```r
ElbowPlot(pbmc)
```

![](https://i.loli.net/2019/12/27/Hp2oFmNcCB3Y4MJ.png)

大多数软件都是通过拾取拐点处的 pc 作为选定数目

### 细胞聚类

```r
pbmc <- FindNeighbors(object = pbmc, dims = 1:10)
pbmc <- FindClusters(object = pbmc, resolution = 0.5)
```

这里的 dims 为上一步计算所用的维度数，而 resolution 参数控制聚类的数目，针对 3K 的细胞数目，最好的范围是**0.4-1.2**

```r
head(Idents(pbmc), 5)
```

```bash
## AAACATACAACCAC AAACATTGAGCTAC AAACATTGATCAGC AAACCGTGCTTCCG AAACCGTGTATGCG
##              1              3              1              2              6
## Levels: 0 1 2 3 4 5 6 7 8
```

### 执行非线性降维

这里注意，这一步聚类有两种聚类方法(umap/tSNE)，两种方法都可以使用，但不要混用，这样，后面的结算结果会将先前的聚类覆盖掉，只能保留一个
本文采用基于 umap 的聚类方法

```r
pbmc <- RunUMAP(object = pbmc, dims = 1:10)
DimPlot(object = pbmc, reduction = "umap")
```

![](https://i.loli.net/2019/12/27/FRckxesLIfd7Mth.png)
完成聚类后，一定要记住保存数据，不然重新计算可要头疼了

```r
saveRDS(pbmc, file = "../output/pbmc_tutorial.rds")
```

### 寻找每个聚类中显著表达的基因

```r
cluster1.markers <- FindMarkers(object = pbmc, ident.1 = 1, min.pct = 0.25)
head(x = cluster1.markers, n = 5)
```

这样是寻找单个聚类中的显著基因

```r
cluster5.markers <- FindMarkers(object = pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(x = cluster5.markers, n = 5)
```

这样寻找所有聚类中显著基因，计算速度很慢，需要等待

另外，我们有多种方法统计基因的显著性

```r
FeaturePlot(object = pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ",
    "PPBP", "CD8A"))
```

![](https://i.loli.net/2019/12/27/ciNovYPm1EIZkVu.png)

```r
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(object = pbmc, features = top10$gene) + NoLegend()
```

![](https://i.loli.net/2019/12/27/xsQuniZyFP7etoV.png)

剩下的便是寻找基因 marker 并对细胞类型进行注释

> 你可能想要知道如何对多样本进行整合分析，请参考[教程：整合刺激性和对照性 PBMC 数据集，以学习细胞类型特异性反应](https://www.hongguangblog.cn/archives/79/)

> 你可能想了解如何灵活操作 Seurat 的 S4 对象，以便轻松的提取表达矩阵；亦或者想要做一些不一样的可视化效果，例如采用平均表达量做热图展示等，请参考[Seurat3.1 的灵活操作指南](https://www.hongguangblog.cn/archives/77/)

> 你可能想对单细胞数据做进一步分析，例如功能分析，请参考[10X 单细胞数据针对细胞及其亚型的基因集功能分析和转录调控分析](https://www.hongguangblog.cn/archives/13/)

## 全自动细胞类型注释

众所周知，细胞类型的注释是最困难的一步，除非你有很强的对细胞基因的敏感度，不然很难识别细胞类型。

通常识别细胞类型的方法主要有三种

- 根据 Marker 基因，采用[CellMarker](http://biocc.hrbmu.edu.cn/CellMarker/)或者[panglaoDB](https://panglaodb.se/)数据库，进行细胞注释,可以采用超几何分布算法来进行精确性验证，方法参考[ClusterProfiler:真的不只是富集分析](https://www.hongguangblog.cn/archives/54/)
- 从文献中获取已经验证的 Marker
- 采用一些自动化注释的软件，例如[`SingleR`](http://www.bioconductor.org/packages/release/bioc/html/SingleR.html),[`Garnett`](https://cole-trapnell-lab.github.io/garnett/docs/#1b-train-your-own-classifier),[`celaref`](http://www.bioconductor.org/packages/release/bioc/html/celaref.html)等

这里简单介绍下`SingleR`

`SingleR`:一个全自动细胞注释的 R 包，用法很简单

### 软件安装

```bash
BiocManager::install("SingleR")
browseVignettes("SingleR")
```

### 创建 SingleR 对象

从头预测的方法请参考[官方教程](http://www.bioconductor.org/packages/release/bioc/vignettes/SingleR/inst/doc/SingleR.html)

因为我们刚刚从 Seurat 过来的，所以我们应该很想知道 Seurat cluster 的细胞注释结果，因此，对 Seurat 的结果进行注释

我们这里采用两个人类的参考集去做细胞注释

```r
library(Seurat)
library(SingleR)
library(dplyr)
library(tibble)
hpca.se <- HumanPrimaryCellAtlasData()
bpe.se <- BlueprintEncodeData()
```

读入`Seurat`对象转换为`SingleCell`支持的对象

```r
seurat.obj <- readRDS("../output/pbmc_tutorial.rds")
seurat.obj@meta.data$cell.type <- Idents(seurat.obj)
test <- as.SingleCellExperiment(seurat.obj)
```

采用两个参考集一起进行注释，

```r
Anno <- SingleR(test = test,
            ref = list(HP = hpca.se , BP = bpe.se),
            labels = list(hpca.se$label.main , bpe.se$label.main),
            method = "cluster",
            cluster = test$cell.type)
```

提取需要的细胞分类信息

```r
Anno$cluster <- rownames(Anno)
fin <- Anno %>% dplyr::tbl_df() %>% dplyr::select(cluster,labels)
```

你也可以将细胞注释信息重新添加到`Seurat`对象中去

```r
new.cluster.ids <- fin$labels
names(new.cluster.ids) <- levels(seurat.obj)
seurat.obj <- RenameIdents(seurat.obj, new.cluster.ids)
```

## 伪时间分析

伪时间分析建议采用 monocle3.0 软件

### 软件安装

```bash
##安装依赖
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor'))
##安装monocle3
devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3')
library(monocle3)
```

### 标准分析流

```r
library(Seurat)
library(monocle3)
endo<-readRDS("../output/pbmc_tutorial.rds")
data <- endo@assays$RNA@counts
pd <-  endo@meta.data
```

不要直接把`meta.data`放入`pd`中去，太多没用的信息了，先清理下

```r
new_pd<-select(pd,nCount_RNA,nFeature_RNA,percent.mt)
new_pd$Cell.type<-Idents(endo)
head(new_pd)
                      nCount_RNA nFeature_RNA percent.mt
AAACCTGGTATCAGTC-1_1      4835         1383   3.536711
AACCGCGCATTCCTGC-1_1      4029         1230   2.655746
AACGTTGAGCCCAATT-1_1      2576          918   3.377329
AAGACCTGTGGTTTCA-1_1      2841         1044   5.455825
AAGTCTGCATGAAGTA-1_1      3731         1213   3.886358
ACAGCCGCATCGGACC-1_1      5915         1918   3.043111
                                     Cell.type
AAACCTGGTATCAGTC-1_1                       DC2
AACCGCGCATTCCTGC-1_1                      MDSC
AACGTTGAGCCCAATT-1_1                      MDSC
AAGACCTGTGGTTTCA-1_1 Dendritic.cells.activated
AAGTCTGCATGAAGTA-1_1                      MDSC
ACAGCCGCATCGGACC-1_1 Dendritic.cells.activated

fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
```

生成`monocle`的`cds`对象

```r
cds <- new_cell_data_set(data,cell_metadata  = new_pd,gene_metadata  = fData)
```

运行下游标准流程，无论如何都得跑，是必须的默认流程

```r
cds <- preprocess_cds(cds, num_dim = 30)
#umap
cds <- reduce_dimension(cds,umap.n_neighbors = 20L)
#color by seurat cluster
plot_cells(cds,label_groups_by_cluster=FALSE,color_cells_by = "Cell.type")
```

![](https://i.loli.net/2019/11/07/GZzaRVeqIoCuEKn.png)

```r
#cluster
cds <- cluster_cells(cds,resolution = 0.5)
#color by monocle cluster
plot_cells(cds, color_cells_by = "partition",label_groups_by_cluster=FALSE)
```

![](https://i.loli.net/2019/11/07/dpSGentaN4j8OvB.png)

```r
cds <- learn_graph(cds)
plot_cells(cds,color_cells_by = "Cell.type",label_groups_by_cluster=FALSE,label_leaves=TRUE,label_branch_points=TRUE)
```

![](https://i.loli.net/2019/11/07/7UcV4zpCFmYLhoB.png)

### 定义时间开始节点&&伪时间分析

```r
##一个有用的寻找起源节点的函数
get_earliest_principal_node <- function(cds, time_bin="Dendritic.cells.activated"){
  cell_ids <- which(colData(cds)[, "Cell.type"] == time_bin)
  closest_vertex <-
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}

#order cell
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
#pseudotime analysis
plot_cells(cds,color_cells_by = "pseudotime",label_cell_groups=FALSE,label_leaves=FALSE,label_branch_points=FALSE,graph_label_size=1.5)
```

![](https://i.loli.net/2019/11/07/fbhs7NUTH4gQZOk.png)

> PS:对多样本进行伪时间分析，请参考[使用 Monocle3 对多样本单细胞数据进行伪时间分析](https://www.hongguangblog.cn/archives/48/)

本文纯属原创，部分数据采用官方教程，转载需标明出处

> 教程链接：
>
> - [教程：整合刺激性和对照性 PBMC 数据集，以学习细胞类型特异性反应](https://www.hongguangblog.cn/archives/79/)
> - [Seurat3.1 的灵活操作指南](https://www.hongguangblog.cn/archives/77/)
> - [10X 单细胞数据针对细胞及其亚型的基因集功能分析和转录调控分析](https://www.hongguangblog.cn/archives/13/)
> - [ClusterProfiler:真的不只是富集分析](https://www.hongguangblog.cn/archives/54/)
> - [使用 Monocle3 对多样本单细胞数据进行伪时间分析](https://www.hongguangblog.cn/archives/48/)
