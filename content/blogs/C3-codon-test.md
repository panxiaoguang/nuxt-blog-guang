---
title: Python编写拆分Barcode的脚本，并用Codon编译为Native code
date: 16th Mar 2023
description: Python编写拆分Barcode的脚本，并用Codon编译为Native code.
image: /blogs-img/codon.webp
alt: Python编写拆分Barcode的脚本，并用Codon编译为Native code
ogImage: /blog-img/codon.webp
tags: ['python']
published: true
---

### 摘要

成对的reads中，read_2的开头包含两份barcode序列，分别长10bp,中间有一段固定长度为15bp的序列分割，例如
`ATCTATGACATGTTACGTTAACTCCNATCTATCACTTAGCGCTGNCCCTGTCCTCTACACTCCACCCCCTCCCCACCAGACTAAACAACGCCCTTTCCCC`

该序列中`ATTTATGACA`及`AATCTATCAA`为barcode序列。要注意，barcode因为测序的原因存在一定的错配，需要对其有一定的容纳。

### 具体过程

#### 1. 编写识别Barcode容错的程序

Barcode和序列本质上是字符串，而最简单的容错就是使用汉明距离来计算出两个等长字符串的错配个数，并加以限制。

```python
def hamming_distance(s1: str, s2: str) -> int:
    """Return the Hamming distance between equal-length sequences"""
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))
```
根据标题所示，我们需要将Python代码编译为本地代码，因此需要对一些变量加以标注，以生成高性能代码。
Codon中对变量类型的声明`s1:str`，在最新版本的python中是存在的。

#### 2. 编写一个解析压缩格式fastq文件的程序

这里因为测序数据量很大，不建议直接读入内存，我们写一个生成器。另外，因为Codon不支持`BioPython`包，因此需要我们自己来写。

```python
import gzip
@tuple
class FastqRecord:
    id: str
    sequence: str
    quality: str

def parse_gzip_fastq(filename:str):
    with gzip.open(filename, 'rt') as f:
        while True:
            try:
                lines = [next(f).strip("\n") for _ in range(4)]
            except:
                break
            if not all(lines):
                break
            id = lines[0][1:]
            sequence = lines[1]
            quality = lines[3]
            yield FastqRecord(id, sequence, quality)
```
#### 3. 读取Barcode序列

Barcode是任意两个序列的2V2组合，而且在read_2上，所以需要反向互补。

```python
import itertools

def reverse_complement(seq: str) -> str:
    """Compute reverse complement of a sequence"""
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    letters = list(seq)
    letters = [basecomplement[base] for base in letters]
    return ''.join(letters[::-1])

def read_barcode(barcode: str):
    """Read barcode from file"""
    codes = List[str]()
    with open(barcode, 'r') as f:
        codes = [reverse_complement(line.strip("\n")) for line in f]
    return [cc for cc in itertools.product(codes, codes)]
```
#### 最终拆分Barcode的Main函数如下

```python
def split_barcodes(fq1:str,fq2:str, barcodes: List[Tuple[str,str]]):
    """make files for every barcode"""
    read1_files = {barcode:open(f"{barcode[0]}_{barcode[1]}_1.fastq","a+") for barcode in barcodes}
    read2_files = {barcode:open(f"{barcode[0]}_{barcode[1]}_2.fastq","a+") for barcode in barcodes}
    """read raw fastq files"""
    for record1, record2 in zip(parse_gzip_fastq(fq1), parse_gzip_fastq(fq2)):
        for barcode in barcodes:
            barcode_l, barcode_r = barcode
            if hamming_distance(record2.sequence[:len(barcode_l)], barcode_l) <= 1 and hamming_distance(record2.sequence[len(barcode_r)+15:len(barcode_r)*2+15], barcode_r) <= 3:
                            read1_files[barcode].write(f"{record1.id}\n")
                            read1_files[barcode].write(record1.sequence+"\n")
                            read1_files[barcode].write("+\n")
                            read1_files[barcode].write(record1.quality+"\n")
                            read2_files[barcode].write(f"{record2.id}\n")
                            read2_files[barcode].write(record2.sequence+"\n")
                            read2_files[barcode].write("+\n")
                            read2_files[barcode].write(record2.quality+"\n")
    [f.close() for f in read1_files.values()]
    [f.close() for f in read2_files.values()]
```
### 编译成本地代码

下载codon编译器
`/bin/bash -c "$(curl -fsSL https://exaloop.io/install.sh)"`

然后执行编译
`codon build -release -o split_barcode split_barcodes.codon`

### 总结

`Codon`作为一款高性能的兼容Python语法的编译器，可以在一定程度上提升`Python`的执行速度到`C/C++`的水平，原则上只要会写`Python`，那么就掌握了99%的`Codon`了，但是在整个过程中我们也发现了，目前`Codon`只支持部分的Python标准库，有些库函数还不健全，例如` GZfile `的`readline()`函数就不存在，需要我们用`next(iter)`来替代，并添加额外的`EOF` 检测。在读取`fastq`文件时，我们也无法调用第三方库，需要自行“造轮子”。虽然提升了执行速度，但是编写代码的速度依旧降低了。

但是，不要忘了现在我们有强大的`chatGPT`了，针对一些功能，你完全可以让`chatGPT`帮你写好代码，这样重复造轮子的工作便可以由AI 取代，这在一定程度上算是减少了`Codon`的缺陷了。