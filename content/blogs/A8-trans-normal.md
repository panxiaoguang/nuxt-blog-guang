---
title: 再说转录组数据标准化（TPM，RPKM，FPKM）
date: 7th Oct 2020
description: 再说转录组数据标准化（TPM，RPKM，FPKM）.
image: /blogs-img/RNA.webp
alt: 再说转录组数据标准化（TPM，RPKM，FPKM）
ogImage: /blog-img/RNA.webp
tags: ['转录组', '生物信息', 'r语言']
published: true
---

### 基础概念讲解
在RNA-Seq的分析中，我们常用**RPKM、FPKM和TPM**作为转录组数据定量的表示方法。

它们都是对表达量进行标准化的方法，为何不直接用read数表示，而选标准化呢?

> 因为落在一个基因区域内的read数目取决于基因长度和测序深度。基因越长read数目越多，测序深度越高,则一个基因对应的read数目也相对越多。所以必须要标准化，而标准化的对象就是基因长度与测序深度。

### RPKM:

Reads Per Kilobase of exon model per Million mapped reads

(每千个碱基的转录每百万映射读取的reads)，主要用来对单端测序（single-end RNA-seq）进行定量的方法。

RPKM(推荐软件，Range) 的计算公式：

$$RPKM=\frac{total\ exon\ reads}{(mapped\ reads\ (Millions)\ \times\ exon\ length(Kb))}$$

**total exon reads**：某个样本mapping到特定基因的外显子上的所有的reads；

**mapped reads ( Millions )** :某个样本的所有reads总和；

**exon length( KB )**：某个基因的长度（外显子的长度的总和，以KB为单位）。

你可以用这个公式计算基因，外显子，转录本的表达，这里以基因的表达为例进行说明。在一个样本中一个基因的RPKM等于落在这个基因上的总的read数(total exon reads)与这个样本的总read数(mapped reads (Millions))和基因长度(exon length( KB )) 的乘积的比值。

### FPKM:

Fragments Per Kilobase of exon model per Million mapped fragments

(每千个碱基的转录每百万映射读取的fragments)，主要是针对pair-end测序表达量进行计算。

FPKM (推荐软件，cufflinks) 和RPKM 的计算方法基本一致。

FPKM和RPKM的区别就是一个是*fragment*，一个是*read*。

对于单末端测序数据，由于Cufflinks计算的时候是将一个read当做一个fragment来算的，故而FPKM等同于RPKM。

对于双末端测序而言，如果一对paired-read都比对上了，那么这一对paired-read称之为一个fragment，而如果一对paired-Read中只有一个比对上了，另外一个没有比对上，那么就将这个比对上的read称之为一个fragment.而计算RPKM时，如果一对paired-read都比对上了会当成两个read计算，而如果一对paired-read中只有一个比对上了，另外一个没有比对上，那么就计read数为1。 故而即使是理论上将各个参数都设置成一样的，也并不能说FPKM=2RPKM。对于单末端测序，虽然理论上FPKM等同于RPKM, 但是实际上即使是使用同一个mapping软件得到的mapping结果，然后再分别去计算同一个基因的RPKM (自己人工计算，或者用现成的一些软件都能算)和FPKM(用Cufflinks计算)，结果却仍然是不同，因为Cufflinks有自己的模型和自己的一些内在算法。

### RPM/CPM:

Reads/Counts of exon model per Million mapped reads (每百万映射读取的reads).

RPM的计算公式：

$$RPM=\ \frac{total\ exon\ reads}{\ mapped\ reads\ (Millions)}$$

**total exon reads**：某个样本mapping到特定基因的外显子上的所有的reads；

**mapped reads (Millions)** :某个样本的所有reads总和；

这就是个占比统计，忽视了转录本长度的影响

### TPM：

Transcripts Per Kilobase of exonmodel per Million mapped reads (每千个碱基的转录每百万映射读取的Transcripts)，优化的RPKM计算方法，可以用于同一物种不同组织的比较。

TPM (推荐软件，RSEM) 的计算公式：

$$TPMi=\;\frac{(\;Ni/Li\;)\ast1000000 {\;SUM(\;Ni/Li+\dots\dots..+\;Nm/Lm\;)}$$


**Ni**：mapping到基因i上的read数；

**Li**：基因i的外显子长度的总和。

在一个样本中一个基因的TPM：先对每个基因的read数用基因的长度进行校正，之后再用校正后的这个基因read数(Ni/Li)与校正后的这个样本的所有read数（sum(Ni/Li+……..+ Nm/Lm)）求商。由此可知，TPM概括了基因的长度、表达量和基因数目。TPM可以用于同一物种不同组织间的比较，因为sum值总是唯一的。


### 该选择哪个作为我的标准化方法？

### TPM一定是万金油？
从概念上我们可以知道，RPKM和FPKM可以优化样本内的基因比较，但是在样本之间比较时，会存在很大的偏见。这个显而易见，因为他们在**不同样本之间的总和都不一致**。

因此，为了可以在多个样本之间比较基因表达量的差异，TPM应运而出，但是TPM真就是万能的吗？

我们可以举一个简单的例子，假设有4个基因，长度都是3bp,readcounts也都是3个，那么，每个基因的TPM=1/4.

假设另一个样本也有这4个基因，但是最后一个基因由于表达量升高而变成了15，那么，前三个基因的TPM=1/8，而最后一个差异基因的TPM=5/8.

这里就引入了一个**偏见**，因为TPM值是相对表达量，本质上仍属于比例的一种，那么这个**相对表达会在总值（分母）改变时发生变化**，那么，我们再比较基因的时候，发现这四个基因都发生了表达量的改变，而实际上只有最后一个是差异基因。

> 那为什么很多文献说TPM可以用来比较多个样本之间的差异呢？我认为他们是默认了样本之间不存在批次效应，即测序深度是一致的或者说相近的，这样保证了分母不会有巨大的变化，那么相对表达就是在样本之间比较差异的金标准。

但是，如果样本之间存在测序深度的差异，那么使用TPM比较样本差异必然会引入偏见，这也是为什么`Deseq2`等差异分析软件不会去选择TPM作为输入了。

### 样本间的标准化方法

假设组间差异成分很大，或者存在很大的批次效应。我们不能用TPM去比较差异，而我们的差异分析软件就会假设我们的**Raw Read Counts**符合负二项分布的模式，假设大部分基因都不是差异基因，差异基因发生在少数基因身上，因此，只有readcount才能符合假设，我们才能用组间矫正的方法例如`edger`包的`TMM`方法。

一个理想的方法便是对TPM执行组间矫正，但是TPM本身并不符合负二项分布的模式，因此会引入偏见。

### TPM如何计算？

* Resem
* Salmon

FeatureCounts的定量结果如何计算TPM？

答案也很简单，我们可以按照公式计算即可，例如这样：

```r
library(dplyr)

df=tibble(gene=c("A","B","C","D"),length=c(12,21,33,45),
          readCount=c(120,13,26,390))

df<-df%>%mutate(Ratio=readCount/length)%>%
     mutate(Sum=sum(Ratio),TPM=Ratio/Sum*1e6)
```
这样就可以简单的计算TPM了，但是这里的并非准确的TPM值，因为我们需要的TPM的length为有效长度而非转录本的长度，而

\EffLength=\feature Length - \average fragment length +1

当然这里的有效长度依然是估计，所以我们需要使用额外的软件去计算插入长度，从而求得片段长度。例如`picard`

最后，我们可以从比对后的BAM文件直接获取基因的定量TPM，通过TPMCalculator软件

文献：[TPMCalculator: One-Step Software to Quantify mRNA Abundance of Genomic Features](https://pubmed.ncbi.nlm.nih.gov/30379987/)
脚本：[Github:TPMCalculator](https://github.com/ncbi/TPMCalculator)