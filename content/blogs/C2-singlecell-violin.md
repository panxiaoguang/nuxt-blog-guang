---
title: 单细胞数据如何绘制stacked violin？
date: 15th Feb 2023
description: 单细胞数据如何绘制stacked violin？
image: /blogs-img/R.webp
alt: 单细胞数据如何绘制stacked violin？
ogImage: /blog-img/R.webp
tags: ['r语言', '单细胞测序']
published: true
---

Python的`Scanpy`包和`Seurat`包一样，是单细胞数据处理的利器，其中，`Scanpy`中有一种堆积的小提琴图，可以很好的展示marker的表达情况，但是在`Seurat`中并没有内置命令。因此，我自己尝试提取数据并用`ggplot2`包来画该图。

首先来展示以下画图的成果，如图

![](https://s3.ax1x.com/2021/01/21/s4ZmSe.png)

那么直接上命令吧！

```r
###载入需要的R包
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
##load测试数据
data<-readRDS("cluster_1.afterclu.rds")
```

```r
##第一个函数从Seurat对象中获取表达值并转化为需要的格式tidy
gotData<-function(seurat.obj,features,groups){
  mat<-GetAssayData(data,assay = "RNA",slot = "data")[features,]
  plotData<-mat%>%
    as.data.frame()%>%
    tibble::rownames_to_column(var="gene")%>%
    as_tibble()%>%
    tidyr::pivot_longer(names_to = "cell",values_to="exp",cols=2:(ncol(mat)+1))
  cellmeta<-data@meta.data%>%
    tibble::rownames_to_column(var="cell")%>%
    as_tibble()%>%
    select(cell,sym(groups))
  plotData<-plotData%>%
    left_join(cellmeta,by="cell")%>%
    setNames(c("gene","cell","value","cellID"))
  plotData
}
##第二个函数画图
plot_stacked_violin<-function(plotData,xlab,ylab,cols){
  ggplot(plotData,aes(y=cellID,x=value,fill=cellID))+
    geom_violin()+
    facet_wrap(.~gene)+
    theme_test()+
    xlab(xlab)+
    ylab(ylab)+
    scale_fill_manual(values = cols)+
    theme(panel.spacing=unit(0,"cm"),
          strip.background = element_rect(fill="transparent",color = "white"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.border = element_rect(size=0.7,colour = "black"),
          strip.text = element_text(size=10,face = "italic"),
          axis.text.y = element_text(size = 11.5,face="bold"),
          axis.title.y = element_text(size = 13))+
    NoLegend()
}
```

下面是使用方法

```r
plot_stacked_violin(gotData(data,c("CD7","KLRB1","NKG7"),"integrated_snn_res.0.1"),"","Cell cluster",c("#8dd3c7","#ffffb3","#bebada","#fb8072"))

```