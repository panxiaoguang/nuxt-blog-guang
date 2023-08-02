---
title: 不同的语言处理gzip压缩文件的时间对比
date: 27th Feb 2023
description: 不同的语言处理gzip压缩文件的时间对比.
image: /blogs-img/language.webp
alt: 不同的语言处理gzip压缩文件的时间对比
ogImage: /blog-img/language.webp
tags: ['r语言', 'python', 'julia']
published: true
---

### 首先在shell中测试如下命令
```bash
#!/bin/sh
time gzip -d -c risearch_chr1:143971112-143971134:+:FAM72C.out.gz > risearch_chr1:143971112-143971134:+:FAM72C.out
```

```
0.28s user 0.02s system 99% cpu 0.297 total
```
### 然后测试python

```python
import gzip
def parse_gzip_py(ris_file):
    inf = gzip.open(ris_file, 'rt')
    with open("test_py_gz.txt","w") as f:
        for line in inf:
            f.write(line+"\n")

    inf.close()
```

```
%timeit parse_gzip_py("risearch_chr1:143971112-143971134:+:FAM72C.out.gz")
399 ms +/- 1.91 ms per loop
```

### Julia 里面有两个解析Gzip的包，分别是`GZip.jl` 和   `CodecZlib.jl`。

我们分别来测试一下

```r
using GZip
function parse_gzip_file(filename::String)
    out = open("test_gzipjl.txt", "w")
    zips = GZip.open(filename)
    while !eof(zips)
        line = readline(zips)
        println(out, line)
    end
    close(out)
end
```

```r
@time parse_gzip_file("risearch_chr1:143971112-143971134:+:FAM72C.out.gz")
0.388812 seconds
```

### 然后是`codeczlib`,必须用另一个包调用它

```r
using TranscodingStreams
using CodecZlib

function parse_gzi_trans(filename::String)
    out = open("test_trans.txt", "w")
    stream = GzipDecompressorStream(open(filename))
    for line in eachline(stream)
        println(out, line)
    end
    close(out)
end
```

```r
@time parse_gzi_trans("risearch_chr1:143971112-143971134:+:FAM72C.out.gz")
0.280360 seconds
```

### 结论：

Julia语言使用GZip包的时候，速度要慢于shell，快于python；
使用CodecZlib的时候，速度快于shell 和 python。但是整体时间来看，最快的比最慢的也就快0.1s，这也就意味着，即使是要解压10000个文件，Julia也就比python快16分钟而已。这个在巨大的解压用时面前，并不算什么。


数据文件可从我的github 获取 
(link)[https://github.com/panxiaoguang/crisproff_jl/tree/main/test_Gzip_parse_time]