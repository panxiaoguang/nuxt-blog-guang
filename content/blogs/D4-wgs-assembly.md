---
title: 使用二代和三代WGS数据进行基因组组装
date: 23th Aug 2023
description: 本教程中，我们的目标是使用第二代和第三代全基因组测序 (WGS) 数据组装细菌基因组。我们将以此为例来探讨WGS数据分析，并探讨测序技术之间的差异。
image: https://picshack.net/ib/EmZbTJiiSx.png
alt: WGS assembly
ogImage: https://picshack.net/ib/EmZbTJiiSx.png
tags: ['生物信息']
published: true
---

本教程中，我们的目标是使用第二代和第三代全基因组测序 (WGS) 数据组装细菌基因组。我们将以此为例来探讨WGS数据分析，并探讨测序技术之间的差异。

### 软件安装

#### Docker镜像

您可以使用以下命令从[DockerHub](https://hub.docker.com/r/yanhui09/mac2023_extra)中提取镜像：

```bash
docker pull yanhui09/mac2023_extra
```

#### 用`mamba`安装

建议使用`mamba`在**独立** `conda`环境中安装软件。

假设您已经在系统中安装了`mamba`软件，您可以使用`mamba`创建一个新`conda`环境：

转到下载的数据目录。_确保您知道自己当前的工作目录位置。_

例如，下载的数据目录位于`/home/username/MAC2023-extera`

```bash
cd /home/username/MAC2023-extera
pwd
```

您将看到下载的数据目录的路径，如下所示。

```bash
/home/username/MAC2023-extera
```

现在您可以使用以下命令创建新`conda`环境并安装软件：

```bash
mamba env create -n wgs1 -f envs/env1.yaml
```

激活环境以进行以下分析。

```bash
conda activate wgs1
```

### 使用演示数据进行探索

WGS 数据可以像以前一样从[MAC2023-extra](https://github.com/yanhui09/MAC2023-extra)获取。

#### 首先用`seqkit`来查看数据

```bash
seqkit stat data/wgs/*.fastq.gz data/wgs/ncbi_pacbio_TL110.fasta
```

预期输出：

```bash
file                              format  type  num_seqs     sum_len    min_len    avg_len    max_len
data/wgs/NXT20x_R1.fastq.gz       FASTQ   DNA    200,000  30,200,000        151        151        151
data/wgs/NXT20x_R2.fastq.gz       FASTQ   DNA    200,000  30,200,000        151        151        151
data/wgs/ont_r10_20x.fastq.gz     FASTQ   DNA      7,862  51,201,670        129    6,512.6     87,688
data/wgs/ncbi_pacbio_TL110.fasta  FASTA   DNA          1   2,566,312  2,566,312  2,566,312  2,566,312
```

这是来自费氏丙酸杆菌菌株的一组测序数据。我们对 Illumina 和 Oxford Nanopore Technologies (ONT) 读数的测序数据进行二次采样，覆盖率达到 20 倍。PacBio 参考基因组来自[NCBI RefSeq 数据库](https://www.ncbi.nlm.nih.gov/nuccore/NZ_CP085641.1)。

> Q1：该菌株的基因组大小是多少？测序覆盖度是如何计算的？
> 
> illumina：30,200,000\*2/2,566,312≈20
> 
> ONT：51,201,670/2,566,312≈20

我们用`fastqc`来看看测序数据的质量。

```bash
mkdir -p fastqc/illumina fastqc/ont_r10
fastqc data/wgs/NXT20x_R*.fastq.gz -o ./fastqc/illumina
fastqc data/wgs/ont_r10_20x.fastq.gz -o ./fastqc/ont_r10 
```

您可以打开`.html`该`fastqc`目录下的文件来查看测序数据的质量。

**ONT**

![](https://picshack.net/ib/4RmGnArPd8.png)

**Illumina**

![](https://picshack.net/ib/WdQeuUuTKy.png)  

> 我们可以看到 ONT 读取比 Illumina 读取更长，但包含更多错误。

#### 使用 Illumina 数据进行基因组组装

##### `trimmomatic`去除接头

`trimmomatic`是一个用于移除接头和低质量数据的工具。[\[阅读更多\]](http://www.usadellab.org/cms/?page=trimmomatic)

使用[Nextera文库制备试剂盒从](https://emea.illumina.com/products/by-type/sequencing-kits/library-prep-kits/nextera-xt-dna.html)[NextSeq](https://emea.illumina.com/systems/sequencing-platforms/nextseq.html)平台收集 Illumina 读数。

```bash
mkdir illumina
trimmomatic PE -threads 4 -phred33 data/wgs/NXT20x_R1.fastq.gz data/wgs/NXT20x_R2.fastq.gz illumina/NXT20x_R1_paired.fastq.gz illumina/NXT20x_R1_unpaired.fastq.gz illumina/NXT20x_R2_paired.fastq.gz illumina/NXT20x_R2_unpaired.fastq.gz ILLUMINACLIP:data/wgs/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 MINLEN:50
```

预期输出：

```bash
illumina/NXT20x_R1_paired.fastq.gz
illumina/NXT20x_R1_unpaired.fastq.gz
illumina/NXT20x_R2_paired.fastq.gz
illumina/NXT20x_R2_unpaired.fastq.gz
```

#### 数据质量控制`bbmap`

`bbmap`是一套用于测序读数质量控制的工具。它可用于删除冗余数据和 _PhiX_ 对照。[\[阅读更多\]](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmap-guide/)

`clumpify.sh`是一个删除冗余数据的工具。[\[阅读更多\]](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/clumpify-guide/)

`bbduk.sh`是一种去除污染数据的工具（例如，宿主基因组、_PhiX_ 对照）。[\[阅读更多\]](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/)

```bash
#clumpify
clumpify.sh in=illumina/NXT20x_R1_paired.fastq.gz in2=illumina/NXT20x_R2_paired.fastq.gz out=illumina/NXT20x_R1_paired_dedup.fastq.gz out2=illumina/NXT20x_R2_paired_dedup.fastq.gz dedupe optical spany adjacent
# bbduk
bbduk.sh in=illumina/NXT20x_R1_paired_dedup.fastq.gz in2=illumina/NXT20x_R2_paired_dedup.fastq.gz out=illumina/NXT20x_R1_paired_dedup_deduk.fastq.gz out2=illumina/NXT20x_R2_paired_dedup_deduk.fastq.gz ref=data/wgs/phiX174.fasta k=31 hdist=1
```

预期输出：

```bash
illumina/NXT20x_R1_paired_dedup.fastq.gz
illumina/NXT20x_R1_paired_dedup_deduk.fastq.gz
illumina/NXT20x_R2_paired_dedup.fastq.gz
illumina/NXT20x_R2_paired_dedup_deduk.fastq.gz
```

##### 基因组组装`spades`

`spades`是一个用于短读长的基因组组装器。[\[阅读更多\]](https://github.com/ablab/spades)

```bash
spades.py --isolate -t 4 -1 illumina/NXT20x_R1_paired_dedup_deduk.fastq.gz -2 illumina/NXT20x_R2_paired_dedup_deduk.fastq.gz -o illumina/spades
```

预计组装：

```bash
illumina/spades/contigs.fasta
```

> _Q2：我们在这里从分离株中组装了细菌基因组。如果我们有宏基因组样本怎么办？_
> 
> 我们可以使用`--meta`选项在`spades`中来组装宏基因组样本。[\[阅读更多\]](https://github.com/ablab/spades#basic-options)

#### 使用 ONT 数据进行基因组组装

##### 使用`guppy`（可选）或者`porechop`移除接头

> 本练习中未提供`guppy`和`porechop`的安装如果您想使用它们，请尝试自行安装。

`guppy`是一种用于 ONT 数据的碱基识别和接头修剪的工具。`guppy`不是开源的，因此需要注册ONT帐户才能查看其文档和下载。[\[阅读更多\]](https://id.customers.nanoporetech.com/app/nanoporetech-customers_myaccount_1/exk2kkmfwpBAaT3WI697/sso/saml?RelayState=https://community.nanoporetech.com/downloads)

`porechop`是一个开源工具，用于 ONT 数据的接头修剪。[\[阅读更多\]](https://github.com/rrwick/Porechop)

默认情况下，条形码将在`guppy`多路分解步骤中被修剪。我们不会对我们的数据重复进行条形码修剪。

如果您想修剪条形码，可以使用以下命令`guppy`。

```bash
guppy_barcoder -i data/wgs/ont_r10_20x.fastq.gz -s ont_r10/ont_r10_20x_barcoded.fastq.gz --barcode_kits EXP-NBD104 --trim_barcodes
```

以及以下命令`porechop`。

```bash
porechop -i data/wgs/ont_r10_20x.fastq.gz -o ont_r10/ont_r10_20x_porechop.fastq.gz --threads 4
```

##### `seqkit`读取质量控制

`seqkit`是一种用于操作测序数据的工具。[\[阅读更多\]](https://bioinf.shenwei.me/seqkit/)这里我们用来`seqkit`去除短读和低质量数据。

```bash
mkdir ont_r10
seqkit seq -j 4 -Q 10 -m 2000 -i data/wgs/ont_r10_20x.fastq.gz -o ont_r10/ont_r10_20x_f.fastq.gz
```

使用`seqkit stat`来检查质量控制之前和之后的 ONT 读数。

```bash
seqkit stat data/wgs/ont_r10_20x.fastq.gz ont_r10/ont_r10_20x_f.fastq.gz 
```

预期输出：

```bash
file                            format  type  num_seqs     sum_len  min_len  avg_len  max_len
data/wgs/ont_r10_20x.fastq.gz   FASTQ   DNA      7,862  51,201,670      129  6,512.6   87,688
ont_r10/ont_r10_20x_f.fastq.gz  FASTQ   DNA      6,561  48,742,515    2,001  7,429.1   87,688 
```

##### `flye`基因组组装

`flye`是 ONT 推荐的长读基因组组装器。[\[阅读更多\]](https://github.com/fenderglass/Flye)

```bash
flye --nano-raw ont_r10/ont_r10_20x_f.fastq.gz --out-dir ont_r10/flye --threads 4 
```

预期的组装输出：

```bash
ont_r10/flye/assembly.fasta 
```

> _Q3：默认_`flye`用于组装基因组。如果我们有宏基因组样本怎么办？
> 
> 我们可以使用`--meta`选项来在`flye`中组装宏基因组样本。[\[阅读更多\]](https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md#-quick-usage)

##### 使用`racon`与`medaka`进行基因组矫正

长读长基因组组装程序通常会产生连续性高但准确性低的基因组草图。需要额外的矫正步骤来提高基因组草图的准确性。但最佳的矫正策略和工具仍在争论中。

`racon`是一种基于图的相似度算法，用于完善长读长基因组组装。[\[阅读更多\]](https://github.com/isovic/racon)。

`medaka`是官方为ONT数据打造的基于神经网络的数据矫正工具。[\[阅读更多\]](https://github.com/nanoporetech/medaka)

最近，ONT 推荐用`medaka`直接矫正`flye`组装后的数据。但结合`racon`和`medaka`仍然是一种常见的做法。

**这里我们以**`medaka`**在例子：**

由于`medaka`是基于神经网络的，选择合适的模型会影响抛光效果。您可以用来`medaka tools list_models`列出所有可用的模型。对于我们的示例，ONT 读数是从[R10.4.1 flowcell](https://store.nanoporetech.com/eu/flow-cell-r10-4-1.html)收集的，并且`guppy`采用`hac`模式来做碱基生成。

```bash
medaka tools list_models
medaka_consensus -i ont_r10/ont_r10_20x_f.fastq.gz -d ont_r10/flye/assembly.fasta -o ont_r10/medaka -t 4 -m r1041_e82_260bps_hac_v4.1.0
```

预期结果：

```bash
ont_r10/medaka/consensus.fasta
```

`racon`和`medaka`混合模式

```bash
mkdir ont_r10/racon
minimap2 -t 4 -x map-ont ont_r10/flye/assembly.fasta ont_r10/ont_r10_20x_f.fastq.gz > ont_r10/racon/flye_assembly.paf
racon -t 4 ont_r10/ont_r10_20x_f.fastq.gz ont_r10/racon/flye_assembly.paf ont_r10/flye/assembly.fasta > ont_r10/racon/racon.fasta
medaka_consensus -i ont_r10/ont_r10_20x_f.fastq.gz -d ont_r10/racon/racon.fasta -o ont_r10/racon/medaka -t 4 -m r1041_e82_260bps_hac_v4.1.0
```

预期的结果：

```bash
ont_r10/racon/racon.fasta
ont_r10/racon/medaka/consensus.fasta
```

#### ONT 和 Illumina 读取的混合组装

混合组装是结合不同测序技术优势的常见做法。这里我们选择两种常用的混合组装策略：

1.  `pilon`：直接使用 ONT 组装作为主干，并用 Illumina 读取对其进行矫正。
    
2.  `unicycler`：短读优先混合组装。
    

##### `pilon`ILLUMINA 读取矫正

`pilon`是一种用短读长优化基因组组装的工具。[\[阅读更多\]](https://github.com/broadinstitute/pilon)

```bash
mkdir -p hybrid/pilon
bwa index ont_r10/medaka/consensus.fasta
bwa mem -t 4 ont_r10/medaka/consensus.fasta illumina/NXT20x_R1_paired_dedup_deduk.fastq.gz illumina/NXT20x_R2_paired_dedup_deduk.fastq.gz | samtools sort -@ 4 -o hybrid/pilon/aligned.bam
samtools index hybrid/pilon/aligned.bam
pilon --genome ont_r10/medaka/consensus.fasta --frags hybrid/pilon/aligned.bam --output hybrid/pilon/pilon --threads 4
```

预期产出：

```bash
hybrid/pilon/pilon.fasta
```

##### 可选：`unicycler`混合组装

`unicycler`可以进行短读优先混合组装。[\[阅读更多\]](https://github.com/rrwick/Unicycler)

**这需要一些时间，因此我们将在示例中跳过这一步。完成其他步骤后，请随时返回此步骤。**

**Linux用户**

```bash
unicycler -l ont_r10/ont_r10_20x_f.fastq.gz -1 illumina/NXT20x_R1_paired_dedup_deduk.fastq.gz -2 illumina/NXT20x_R2_paired_dedup_deduk.fastq.gz -o hybrid/unicycler --threads 4
```

**MacOS 用户**

由于on的最新`MacOS`版本落后于[0.5.0](https://github.com/rrwick/Unicycler/releases/tag/v0.5.0)，我们附加两个标志和来保持最大的兼容性。`unicyclerconda--no_correct--no_pilon`

```bash
unicycler --no_correct --no_pilon -l ont_r10/ont_r10_20x_f.fastq.gz -1 illumina/NXT20x_R1_paired_dedup_deduk.fastq.gz -2 illumina/NXT20x_R2_paired_dedup_deduk.fastq.gz -o hybrid/unicycler
```

预期的产出：

```bash
hybrid/unicycler/assembly.fasta
```

#### 使用`quast对`组装基因组进行质量评估

我们已经生成了许多具有不同策略的程序集。让我们使用基于PacBio深度测序的完整组装作为参考基因组来评估每个数据集的质量。如果您尚未完成所有组装步骤，您可以使用`./data/wgs/assemblies`提供的组装文件。

```bash
seqkit stat data/wgs/assemblies/*.fasta data/wgs/ncbi_pacbio_TL110.fasta
```

预期输出：

```bash
file                                               format  type  num_seqs    sum_len    min_len      avg_len    max_len
data/wgs/assemblies/hybrid_pilon.fasta             FASTA   DNA          2  2,579,926     29,242    1,289,963  2,550,684
data/wgs/assemblies/hybrid_pilon_proovframe.fasta  FASTA   DNA          2  2,579,979     29,245  1,289,989.5  2,550,734
data/wgs/assemblies/hybrid_unicycler.fasta         FASTA   DNA          1  2,564,177  2,564,177    2,564,177  2,564,177
data/wgs/assemblies/illumina.fasta                 FASTA   DNA        217  2,527,918         78     11,649.4     60,978
data/wgs/assemblies/ont_flye.fasta                 FASTA   DNA          2  2,579,333     29,239  1,289,666.5  2,550,094
data/wgs/assemblies/ont_flye_medaka.fasta          FASTA   DNA          2  2,579,052     29,238    1,289,526  2,549,814
data/wgs/assemblies/ont_flye_racon.fasta           FASTA   DNA          2  2,578,932     29,231    1,289,466  2,549,701
data/wgs/assemblies/ont_flye_racon_medaka.fasta    FASTA   DNA          2  2,578,566     29,229    1,289,283  2,549,337
data/wgs/ncbi_pacbio_TL110.fasta                   FASTA   DNA          1  2,566,312  2,566,312    2,566,312  2,566,312
```

> 与参考相比，各个组装的总长度相似。
> 
> 但序列数 ( `num_seqs`) 和最大序列长度 ( `max_len`) 变化很大。

让我们用`quast`来检查各个组装的质量。

```bash
quast data/wgs/assemblies/*.fasta -r data/wgs/ncbi_pacbio_TL110.fasta -o quast
```

您可以打开`quast`目录中的`report.html`文件来查看各个组装数据集的质量。

![](https://picshack.net/ib/d6khde5lGw.png)

> _Q4：根据示例，您认为哪种测序技术在基因组完整性和连续性方面效果更好？_
> 
> 检查`Genome fraction`和`NGA50`
> 
> _Q5：根据示例，您认为哪种装配的精度最好？_
> 
> 检查`Misassemblies`和`Mismatches`
> 
> _Q6：根据示例，ONT 组件中的主要错误类型是什么？_
> 
> 检查`mismatches`和`indels`。[\[阅读更多\]](https://genome.sph.umich.edu/wiki/Indel)
> 
> _Q7：您认为此示例的最佳组装策略是什么？_

#### 可选：参考引导校正`proovframe`

高频率`indels`会导致 ONT 组装的编码序列 (CDS) 发生移码。通过参考引导校正，我们可以校正 ONT 组装的 CDS 中的移码。

**此步骤是可选的，如果您没有时间，可以跳过它。**

为了使用该工具，由于依赖冲突，`proovframe`我们需要创建另一个`conda`环境。

```bash
conda deactivate
mamba env create -n wgs2 -f envs/env2.yaml
conda activate wgs2
```

##### 基因组注释`prokka`

`prokka`是快速注释细菌基因组的常用工具。[\[阅读更多\]](https://github.com/tseemann/prokka)

```bash
mkdir proovframe
conda activate wgs2
prokka --outdir proovframe/prokka --prefix pacbio --cpus 4 data/wgs/ncbi_pacbio_TL110.fasta
```

预期输出：

```bash
proovframe/prokka/pacbio.err
proovframe/prokka/pacbio.faa
proovframe/prokka/pacbio.ffn
proovframe/prokka/pacbio.fna
proovframe/prokka/pacbio.fsa
proovframe/prokka/pacbio.gbk
proovframe/prokka/pacbio.gff
proovframe/prokka/pacbio.log
proovframe/prokka/pacbio.sqn
proovframe/prokka/pacbio.tbl
proovframe/prokka/pacbio.tsv
proovframe/prokka/pacbio.txt
```

我们有许多来自`prokka`的输出文件。这里我们仅使用翻译后的蛋白质序列（`./proovframe/prokka/pacbio.faa`文件）

##### 移码校正`proovframe`

![](https://picshack.net/ib/Z55fiIPp1W.png)

`proovframe`是 CDS 中基于参考序列的移码校正工具。[\[阅读更多\]](https://www.biorxiv.org/content/10.1101/2021.08.23.457338v1)

```bash
proovframe map -a proovframe/prokka/pacbio.faa -o proovframe/pilon.tsv hybrid/pilon/pilon.fasta
proovframe fix -o proovframe/pilon_corrected.fasta hybrid/pilon/pilon.fasta proovframe/pilon.tsv
```

预期的输出：

```bash
proovframe/pilon_corrected.fasta
```

> _Q8：在之前的报告中，_`shybrid_pilon_proovframe.fasta`_中的_`N'`_会增加。你觉得这是好是坏？_

> 本文为学习记录，课程作者为：[YanHui](https://www.yanh.org/)