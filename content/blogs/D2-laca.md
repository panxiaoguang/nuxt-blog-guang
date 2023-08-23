---
title: 使用LACA从长扩增子中从头挑选OTU
date: 23th Aug 2023
description: LACA是用于长扩增子一致性分析（例如 16S rRNA 基因扩增子分析）的可重复且可扩展的工作流程。它使用用snakemake管理工作流程以及conda来管理环境。
image: https://picshack.net/ib/EmZbTJiiSx.png
alt: LACA selects OTUs
ogImage: https://picshack.net/ib/EmZbTJiiSx.png
tags: ['生物信息']
published: true
---

`LACA`是用于长扩增子一致性分析（例如 16S rRNA 基因扩增子分析）的可重复且可扩展的工作流程。它使用用`snakemake`管理工作流程以及`conda来`管理环境。

### LACA的安装

完整的安装指南`LACA`可[在此处](https://github.com/yanhui09/laca)获取。

您可以根据您的喜好选择使用`docker`或从`GitHub`存储库安装`LACA`

#### Docker镜像

最简单的使用方法是从[Docker Hub](https://hub.docker.com/r/yanhui09/laca)拉取`LACA`镜像以获得跨平台支持

```bash
docker pull yanhui09/laca 
```

> `LACA`是通过`docker`为`linux/amd64`平台而构建的，
> 
> `MacOS`用户需要使用docker容器来运行`LACA`。

#### 从 GitHub 存储库安装

**1.克隆Github仓库并创建隔离`conda`环境**

```bash
git clone https://github.com/yanhui09/laca.git
cd laca
mamba env create -n laca -f env.yaml 
```

**2.安装`LACA`**

为了避免不一致，建议在上面建立的`conda`环境中安装`LACA`

```bash
conda activate laca
pip install --editable . 
```

### 使用LACA运行演示数据

[在这里](https://github.com/yanhui09/laca#usage)找到完整的使用指南。

#### 快速启动示例

```bash
laca init -b /path/to/basecalled_fastqs -d /path/to/database    # init config file and check
laca run all                                         # start analysis 
```

#### 熟悉`LACA`使用

`LACA`很容易使用。您可以使用`laca init`和`laca run`分两步开始新的分析。

如果`LACA`安装在`conda`环境中，请记住激活`conda`环境。

```bash
conda activate laca
laca -h 
```

**要使用 docker 镜像，您需要将数据目录（例如`pwd`）挂载到容器中`/home` 目录。**

```bash
docker run -it -v `pwd`:/home --privileged yanhui09/laca
laca -h 
```

**1.初始化配置文件`laca init`**

`laca init`会在工作目录中生成一个配置文件，其中包含运行`LACA`所需的所有参数。

```bash
laca init -h 
```

**2.开始`laca run`分析**

`laca run`将相应地触发完整的工作流程或定义资源下的特定模块。使用`laca run -h`获得试运行概述。

```bash
laca run -h 
```

#### 使用演示数据集运行`LACA`

**0.确保您已从[此处](https://github.com/yanhui09/MAC2023-extra)下载所需的演示数据集。然后`cd`进入目录。**

例如，输入绝对路径（“长路径”）`/home/me/MAC2023-extra`。

```bash
cd /home/me/MAC2023-extra 
```

如果您尚未下载数据用`Git`下载，

```bash
git clone https://github.com/yanhui09/MAC2023-extra.git
cd ./MAC2023-extra 
```

**1.检查您所在的位置并尝试`laca init`检查生成的`config.yaml`文件。**

```bash
pwd
laca init -b ./data/ont16s -d ./database -w ./laca_output --fqs-min 50
cat ./laca_output/config.yaml 
```

**2.`LACA`伪运行和真实运行**

```bash
laca run all -w ./laca_output -n 
laca run kmerCon -j 4 -w ./laca_output 
```

`LACA`能够生成`otu table`，`taxonomy table`以及`phylogenetic tree`如果您使用`laca run all`运行完整的工作流程。但第一次使用需要时间准备数据库和安装。

作为一个例子，这里我们只运行模块`kmerCon`来根据 kmer 频率提取一致序列。

看看这些共有序列，取第一个序列对`rRNA/ITS`数据库进行[BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)搜索。

```bash
head -n2 ./laca_output/kmerCon/kmerCon.fna 
```

预期输出：

```bash
>pooled_0b000_0cand1
CACAATGGGCGCAAGCCTGATGCAGCGACGCCGCGTGCGGGATGACGGCCTTCGGGTTGTAAACCGCTTTTGACTGGGAGCAAGCCCTTCGGGGTGAGTGTACCTTTCGAATAAGCACCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGTGCAAGCGTTATCCGGAATTATTGGGCGTAAAGGGCTCGTAGGCGGTTCGTCGCGTCCGGTGTGAAAGTCCATCGCTTAACGGTGGATCCGCGCCGGGTACGGGCGGGCTTGAGTGCGGTAGGGGAGACTGGAATTCCCGGTGTAACGGTGGAATGTGTAGATATCGGGAAGAACACCAATGGCGAAGGCAGGTCTCTGGGCCGTCACTGACGCTGAGGAGCGAAAGCGTGGGGAGCGAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGGTGGATGCTGGATGTGGGGACCATTCCACGGTCTCCGTGTCGGAGCCAACGCGTTAAGCATCCCGCCTGGGGAGTACGGCCGCAAGGCTAAAACTCAAAGAAATTGACGGGGGCCCGCACAAGCGGCGGAGCATGCGGATTAATTCGATGCAACGCGAAGAACCTTACCTGGGCTTGACATGTTCCCGACAGCCGTAGAGATACGGCCTCCCTTCGGGGCGGGTTCACAGGTGGTGCATGGTCGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTCGCCCTGTGTTGCCAGCACGTCGTGGTGGGAACTCACGGGGGACCGCCGGGGTCAACTCGGAGGAAGGTGGGGATGACGTCAGATCATCATGCCCCTTACGTCCAGGGCTTCACGCATGCTACAATGGCCGGTACAACGGGATGCGACCTCGCGAGGGGGAGCGGATCCCTTAAAACCGGTCTCAGTTCGGATTGGAGTCTGCAACCCGACTCCATGAAGGCGGAGTCGCTAGTAATCGCGGATCAGCAACGCCGCGGTGAATGCGTTCCCGGGCC 
```

![](https://picshack.net/ib/ERENmtCCeN.jpg)

结果`BLAST`表明，该序列是一个`16S rRNA`基因片段，`Bifidobacterium`一性性超过99%。

> 读取一个 ONT 并进行相同的`BLAST`搜索。您希望看到什么？
> 
> ```bash
> zcat ./data/ont16s/*.fastq.gz | head
> ```
> 
> 本文为学习记录，课程作者为：[YanHui](https://www.yanh.org/)