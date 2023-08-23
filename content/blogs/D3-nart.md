---
title: 使用NART通过读取分类进行长扩增子分析
date: 23th Aug 2023
description: NART设计用于基于图谱的纳米孔扩增子（实时）分析，例如 16S rRNA 基因。NART由NART（Nanopore Amplicon Real-Time entry）和 NAWF（Nanopore Amplicon snakemake WorkFlow entry）组成。通过基于映射的策略提供从基础调用读取到最终计数矩阵的（实时）端到端解决方案。
image: https://picshack.net/ib/EmZbTJiiSx.png
alt: NART analysis
ogImage: https://picshack.net/ib/EmZbTJiiSx.png
tags: ['生物信息']
published: true
---

`NART`设计用于基于图谱的纳米孔扩增子（**实时**）分析，例如 16S rRNA 基因。`NART`由`NART`（Nanopore Amplicon Real-Time entry）和 `NAWF`（Nanopore Amplicon `snakemake` WorkFlow entry）组成。通过基于映射的策略提供从基础调用读取到最终计数矩阵的（实时）端到端解决方案。

`nawf`提供三个选项（`emu`，`minimap2lca`和`blast2lca`）来确定微生物组成。

![](https://picshack.net/ib/O8BXIaIwax.png)

### NART安装

完整的安装指南`NART`可[在此处](https://github.com/yanhui09/nart)获取。

您可以根据您的喜好选择`NART`使用`docker`映像安装（这只是`MacOS`用户的解决方案）或从`GitHub`存储库安装。

#### Docker镜像

最简单的使用方法是从[Docker Hub](https://hub.docker.com/r/yanhui09/nart)拉取`NART`镜像以获得跨平台支持。

```bash
docker pull yanhui09/nart 
```

> `NART`是通过`docker`为`linux/amd64`平台而构建的，
> 
> `MacOS`用户需要使用docker容器来运行`NART`。

#### 从 GitHub 存储库安装

**1.克隆Github仓库并创建隔离`conda`环境**

```bash
git clone https://github.com/yanhui09/nart.git
cd nart
mamba env create -n nart -f env.yaml 
```

**2.安装`NART`**

为了避免不一致，建议`NART`在上述`conda`环境中安装

```bash
conda activate nart
pip install --editable . 
```

### 使用NART运行的演示

[在这里](https://github.com/yanhui09/nart#usage)找到完整的使用指南。

[可以在此处](https://www.youtube.com/watch?v=TkdJGLOscPg)找到视频教程。

#### 快速启动示例

##### 单批次扩增子分析

`nawf`可用于分析Nanopore 运行或批次中的任何单个碱基调用文件。

```bash
nawf config -b /path/to/basecall_fastq -d /path/to/database    # init config file and check
nawf run all                                                   # start analysis 
```

##### 实时分析

`nart`提供实用程序来记录、处理和分析连续生成的`fastq`批次。

在开始实时分析之前，您需要`nawf`根据需要配置工作流程。

```bash
nawf config -d /path/to/database                               # init config file and check 
```

在常见情况下，您需要三个独立的会话来分别处理监控、处理和可视化。

**1.监听 bascall 输出并记录**

```bash
nart monitor -q /path/to/basecall_fastq_dir                    # monitor basecall output 
```

**2.开始新的 fastq 的扩增子分析**

```bash
nart run -t 10                                                 # real-time process in batches 
```

**3.更新特征表以在浏览器中交互式可视化**

```bash
nart visual                                                    # interactive visualization 
```

#### 熟悉`NART`使用

`NART`由两组脚本组成：`nart`和`nawf`，分别控制实时分析和工作流性能。

如果`NART`安装在`conda`环境中，请记住激活环境`conda`。

```bash
conda activate nart
nawf -h
nart -h 
```

**要使用 docker 镜像**，您需要将数据目录（例如`pwd` ）挂载到容器中`/home目录`。

```bash
docker run -it -v `pwd`:/home --network host --privileged yanhui09/nart
nawf -h
nart -h 
```

> **注意：**`--network host`需要`nart monitor`才能正常工作。

#### 使用演示数据集运行`NART`

**0.确保您已从[此处](https://github.com/yanhui09/MAC2023-extra)下载所需的演示数据集。然后进入目录`cd`。**

例如，输入绝对路径（“长路径”）的目录是`/home/me/MAC2023-extra`。

```bash
cd /home/me/MAC2023-extra 
```

如果您尚未下载使用`Git`下载数据，

```bash
git clone https://github.com/yanhui09/MAC2023-extra.git
cd ./MAC2023-extra 
```

####  `nawf`分析已完成的 ONT

**1.1. 检查您所在的位置并尝试`laca init`检查生成的`config.yaml`文件。**

```bash
pwd
nawf config -b ./data/ont16s/*.fastq.gz -d ./database -w ./nart_output
cat ./nart_output/config.yaml 
```

**1.2.开始`nawf`试运行** 

```bash
nawf run all -w ./nart_output -n 
```

##### 2. 实时分析`nart`

**2.1.** 重新生成`config.yaml`不带`-b`标志的文件

```bash
rm -rf ./nart_output
nawf config -d ./database -w ./nart_output
head ./nart_output/config.yaml 
```

> `basecall_fq`检查文件中的更改`config.yaml`。

**2.2.** 监控 bascalling 输出并记录

```bash
nart monitor -q ./data/ont16s -w ./nart_output 
```

`nart monitor`在工作目录中创建一个`fqs.txt`来记录即将到来的`fastq`文件。

**2.3.** 开始扩增子分析（需要新终端）

```bash
nart run -t 4 -w ./nart_output 
```

> 在新终端中，检查以下操作`nart monitor`
> 
> ```bash
> ls ./nart_output
> cat ./nart_output/fqs.txt
> cp ./data/ont16s/*.fastq.gz ./data/ont16s/new.fastq.gz
> cat ./nart_output/fqs.txt 
> ```
> 
> 检查文件内容`fqs.txt`的变化。

`nart run`将特定于批次的计数矩阵存储在文件夹`batches`下。当创建新矩阵时，合并表 ( `otu_table.tsv`) 会迭代更新。

当您在输出文件夹中看到一个otu\_table.tsv，例如`nart_output/`，您可以尝试下面的交互式可视化。

**2.4.** 浏览器中的交互式可视化（需要新终端）

```bash
nart visual -i ./nart_output 
```

在浏览器中打开生成的链接。您预计会看到如下所示的交互式条形图。![](https://picshack.net/ib/ScRCwcE8ur.png)

> `MacOS`用户无法通过`docker`体验`nart visual`。:(
> 
> 主机网络驱动程序仅适用于 Linux 主机，在 Docker Desktop for Mac、Docker Desktop for Windows 或 Docker EE for Windows Server 上不受支持。[\[阅读更多\]](https://docs.docker.com/network/drivers/host/)

本文为教程整理，原作者为[YanHui](https://www.yanh.org/)