---
title: 教程：整合刺激性和对照性PBMC数据集，以学习细胞类型特异性反应
date: 20th Nov 2019
description: 教程：整合刺激性和对照性PBMC数据集，以学习细胞类型特异性反应.
image: /blogs-img/singlecell.webp
alt: 教程：整合刺激性和对照性PBMC数据集，以学习细胞类型特异性反应
ogImage: /blogs-img/singlecell.webp
tags: ['单细胞测序', '生物信息', 'r语言']
published: true
---

本教程介绍了[Kang等人（2017)](https://www.nature.com/articles/nbt.4042)的两组PBMC的对齐方式。在该实验中，将PBMC分为刺激组和对照组，并用干扰素β治疗刺激组。对干扰素的反应引起细胞类型特异性基因表达的变化，这使得对所有数据进行联合分析变得困难，并且细胞按刺激条件和细胞类型聚类。在这里，我们证明了我们的整合策略，如[Stuart和Butler等人（2018年）](https://www.biorxiv.org/content/early/2018/11/02/460147)所述，用于执行整合分析以促进常见细胞类型的鉴定并进行比较分析。尽管此示例演示了两个数据集（条件）的集成，但这些方法已扩展到多个数据集。这个[工作流程](https://satijalab.org/seurat/pancreas_integration_label_transfer.html)提供了整合四个胰岛数据集的示例。


### 整合目标

以下教程旨在概述使用Seurat集成过程可能进行的复杂细胞类型的比较分析。在这里，我们解决了三个主要目标：

* 识别两个数据集中都存在的单元格类型
* 获得在对照和刺激细胞中均保守的细胞类型标记
* 比较数据集以找到对刺激的细胞类型特异性反应

### 设置Seurat对象

为方便起见，我们通过`SeuratData`软件包分发此数据集。
```r
library(Seurat)
library(SeuratData)
library(cowplot)
InstallData("ifnb")
data("ifnb")
ifnb.list <- SplitObject(ifnb, split.by = "stim")

ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
```

### 执行整合

然后，我们使用`FindIntegrationAnchors`函数来识别锚点，该函数将`Seurat`对象的列表作为输入，并使用这些锚点将两个数据集与集成在一起。
```r
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, dims = 1:20)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)
```

### 进行综合分析

现在，我们可以在所有单元上运行单个集成分析！
```r
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)
# Visualization
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
```
![](https://i.loli.net/2019/12/11/J7YB4t9FGzNgjUx.png)
为了并排可视化这两个条件，我们可以使用split.by参数来显示每个以聚类着色的条件。
```r
DimPlot(immune.combined, reduction = "umap", split.by = "stim")
```
![](https://i.loli.net/2019/12/11/qKrMcbdkETNV81n.png)

### 识别保守的细胞类型标记

为了确定跨条件保守的规范细胞类型标记基因，我们提供了该`FindConservedMarkers`功能。此功能对每个数据集/组执行差异基因表达测试，并使用`MetaDE`R软件包中的荟萃分析方法组合p值。例如，无论簇7中的刺激条件如何，我们都可以计算出保守标记的基因（NK细胞）。
```r
DefaultAssay(immune.combined) <- "RNA"
nk.markers <- FindConservedMarkers(immune.combined, ident.1 = 7, grouping.var = "stim", verbose = FALSE)
head(nk.markers)
##           CTRL_p_val CTRL_avg_logFC CTRL_pct.1 CTRL_pct.2 CTRL_p_val_adj
## SNHG12 1.059703e-193      1.3678805      0.335      0.019  1.489200e-189
## HSPH1  2.552539e-139      2.0586178      0.553      0.100  3.587083e-135
## NR4A2  8.671555e-136      1.5545769      0.296      0.023  1.218614e-131
## SRSF2  1.556024e-113      1.6410606      0.704      0.220  2.186680e-109
## BATF    1.573042e-09      0.5991204      0.116      0.042   2.210596e-05
## CD69   1.188324e-122      1.8357378      0.525      0.096  1.669952e-118
##           STIM_p_val STIM_avg_logFC STIM_pct.1 STIM_pct.2 STIM_p_val_adj
## SNHG12 8.090842e-159       1.054494      0.256      0.015  1.137006e-154
## HSPH1   4.097380e-89       1.580183      0.471      0.114   5.758049e-85
## NR4A2   3.130700e-78       1.009556      0.172      0.015   4.399572e-74
## SRSF2  1.829674e-128       1.625081      0.675      0.182  2.571241e-124
## BATF   9.234006e-126       1.354443      0.305      0.031  1.297655e-121
## CD69    3.733167e-78       1.677068      0.688      0.291   5.246220e-74
##             max_pval minimump_p_val
## SNHG12 8.090842e-159  2.119406e-193
## HSPH1   4.097380e-89  5.105078e-139
## NR4A2   3.130700e-78  1.734311e-135
## SRSF2  1.556024e-113  3.659348e-128
## BATF    1.573042e-09  1.846801e-125
## CD69    3.733167e-78  2.376649e-122
```

我们可以为每个簇探索这些标记基因，并使用它们将我们的簇注释为特定的细胞类型。
```r
FeaturePlot(immune.combined, features = c("CD3D", "SELL", "CREM", "CD8A", "GNLY", "CD79A", "FCGR3A", 
    "CCL2", "PPBP"), min.cutoff = "q9")
```

![](https://i.loli.net/2019/12/11/LSz7HRaD1AovJQB.png)

```r
immune.combined <- RenameIdents(immune.combined, `0` = "CD14 Mono", `1` = "CD4 Naive T", `2` = "CD4 Memory T", 
    `3` = "CD16 Mono", `4` = "B", `5` = "CD8 T", `6` = "T activated", `7` = "NK", `8` = "DC", `9` = "B Activated", 
    `10` = "Mk", `11` = "pDC", `12` = "Eryth", `13` = "Mono/Mk Doublets")

DimPlot(immune.combined, label = TRUE)
```
![](https://i.loli.net/2019/12/11/N9qipvX6zWoeUaS.png)

`DotPlot`带有`split.by`参数的函数可用于查看各种条件下保守的细胞类型标记，显示表达水平和表达任何给定基因的簇中细胞的百分比。在这里，我们为13个簇中的每个簇绘制2-3个强标记基因。
```r
Idents(immune.combined) <- factor(Idents(immune.combined), levels = c("Mono/Mk Doublets", "pDC", 
    "Eryth", "Mk", "DC", "CD14 Mono", "CD16 Mono", "B Activated", "B", "CD8 T", "NK", "T activated", 
    "CD4 Naive T", "CD4 Memory T"))
markers.to.plot <- c("CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY", "NKG7", "CCL5", 
    "CD8A", "MS4A1", "CD79A", "MIR155HG", "NME1", "FCGR3A", "VMO1", "CCL2", "S100A9", "HLA-DQA1", 
    "GPR183", "PPBP", "GNG11", "HBA2", "HBB", "TSPAN13", "IL3RA", "IGJ")
DotPlot(immune.combined, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 8, 
    split.by = "stim") + RotatedAxis()
```
![](https://i.loli.net/2019/12/11/OhXroNgvKuTDbm7.png)

### 跨条件鉴定差异表达基因
现在我们已经排列了刺激细胞和对照细胞，我们可以开始进行比较分析，并观察刺激引起的差异。广泛观察这些变化的一种方法是绘制受刺激细胞和对照细胞的平均表达，并在散点图上寻找视觉异常值的基因。在这里，我们采用受刺激的和对照的原始T细胞和CD14单核细胞群体的平均表达，并生成散点图，突出显示对干扰素刺激表现出戏剧性反应的基因。
```r
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
t.cells <- subset(immune.combined, idents = "CD4 Naive T")
Idents(t.cells) <- "stim"
avg.t.cells <- log1p(AverageExpression(t.cells, verbose = FALSE)$RNA)
avg.t.cells$gene <- rownames(avg.t.cells)

cd14.mono <- subset(immune.combined, idents = "CD14 Mono")
Idents(cd14.mono) <- "stim"
avg.cd14.mono <- log1p(AverageExpression(cd14.mono, verbose = FALSE)$RNA)
avg.cd14.mono$gene <- rownames(avg.cd14.mono)

genes.to.label = c("ISG15", "LY6E", "IFI6", "ISG20", "MX1", "IFIT2", "IFIT1", "CXCL10", "CCL8")
p1 <- ggplot(avg.t.cells, aes(CTRL, STIM)) + geom_point() + ggtitle("CD4 Naive T Cells")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
p2 <- ggplot(avg.cd14.mono, aes(CTRL, STIM)) + geom_point() + ggtitle("CD14 Monocytes")
p2 <- LabelPoints(plot = p2, points = genes.to.label, repel = TRUE)
plot_grid(p1, p2)
```
![](https://i.loli.net/2019/12/11/bPnJlFC4tY1xcu9.png)

如您所见，许多相同的基因在这两种细胞类型中均被上调，可能代表保守的干扰素应答途径。

因为我们有信心在各种情况下都能识别出常见的细胞类型，所以我们可以询问相同类型的细胞在不同条件下会改变哪些基因。首先，我们在`meta.data`插槽中创建一列，以保存细胞类型和刺激信息，并将当前标识切换到该列。然后，我们使用它`FindMarkers`来找到受刺激的B细胞和对照B细胞之间不同的基因。请注意，此处显示的许多顶级基因与我们之前绘制的核心干扰素应答基因相同。此外，我们看到的诸如CXCL10的基因对单核细胞和B细胞干扰素的反应也具有特异性，在该列表中也显示出很高的意义。
```r
immune.combined$celltype.stim <- paste(Idents(immune.combined), immune.combined$stim, sep = "_")
immune.combined$celltype <- Idents(immune.combined)
Idents(immune.combined) <- "celltype.stim"
b.interferon.response <- FindMarkers(immune.combined, ident.1 = "B_STIM", ident.2 = "B_CTRL", verbose = FALSE)
head(b.interferon.response, n = 15)
##                 p_val avg_logFC pct.1 pct.2     p_val_adj
## ISG15   9.008631e-168 3.2061225 0.998 0.235 1.265983e-163
## IFIT3   1.566158e-161 3.1240285 0.961 0.049 2.200921e-157
## ISG20   1.206176e-158 2.0549983 1.000 0.662 1.695038e-154
## IFI6    1.544804e-158 2.9139826 0.959 0.077 2.170914e-154
## IFIT1   1.531805e-144 2.8529052 0.899 0.031 2.152645e-140
## MX1     3.490304e-129 2.2712525 0.902 0.113 4.904924e-125
## LY6E    3.320706e-127 2.1767293 0.897 0.146 4.666588e-123
## TNFSF10 2.108526e-114 2.5996223 0.773 0.021 2.963111e-110
## IFIT2   7.373988e-113 2.5114519 0.781 0.035 1.036267e-108
## B2M     1.570724e-101 0.4194467 1.000 1.000  2.207338e-97
## PLSCR1  8.204551e-101 1.9635626 0.792 0.115  1.152986e-96
## IRF7    1.517587e-100 1.8203577 0.837 0.181  2.132665e-96
## CXCL10   7.783861e-92 3.6869595 0.660 0.012  1.093866e-87
## UBE2L6   1.659777e-88 1.4951634 0.855 0.296  2.332484e-84
## EPSTI1   1.324448e-82 1.7540298 0.717 0.103  1.861247e-78
```

可视化基因表达中这些变化的另一种有用方法是`split.by`选择`FeaturePlot`或`VlnPlot`功能。这将显示给定基因列表的`FeaturePlots`，并按分组变量（此处为刺激条件）进行划分。诸如CD3D和GNLY之类的基因是典型的细胞类型标记（对于T细胞和NK / CD8 T细胞），实际上不受干扰素刺激的影响，并且在对照组和受刺激组中显示出相似的基因表达模式。另一方面，IFI6和ISG15是核心干扰素应答基因，因此在所有细胞类型中均被上调。最后，CD14和CXCL10是显示细胞类型特异性干扰素应答的基因。CD14单核细胞受刺激后，CD14表达下降，这可能导致在有监督的分析框架中进行错误分类，从而强调了整合分析的价值。
```r
FeaturePlot(immune.combined, features = c("CD3D", "GNLY", "IFI6"), split.by = "stim", max.cutoff = 3, 
    cols = c("grey", "red"))
```
![](https://i.loli.net/2019/12/11/9Co6NOnejimagsc.png)
```r
plots <- VlnPlot(immune.combined, features = c("LYZ", "ISG15", "CXCL10"), split.by = "stim", group.by = "celltype", 
    pt.size = 0, combine = FALSE)
CombinePlots(plots = plots, ncol = 1)
```
![](https://i.loli.net/2019/12/11/DIyLP7VsZEXKUTf.png)
