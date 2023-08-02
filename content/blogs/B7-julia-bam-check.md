---
title: Julia 短小代码批量检测BAM文件的完整性
date: 15th Feb 2021
description: Julia 短小代码批量检测BAM文件的完整性.
image: /blogs-img/julia.webp
alt: Julia 短小代码批量检测BAM文件的完整性
ogImage: /blog-img/julia.webp
tags: ['julia']
published: true
---

我们在运行bwa mem比对的时候，由于某些不明的原因会造成程序中断，例如内存超了，IO错误，计算节点崩溃等，然而BAM是否完整很难察觉，最终导致后续流程无法运行。这里，我们通过一段简短的代码来检查BAM文件的完整性，代码如下：


```r
function check_eofs(fs::String)
    open(fs, "r") do IO
        seekend(IO)
        eof_positions = position(IO) - 28
        seek(IO, eof_positions)
        eof_markers = UInt8[0x1f, 0x8b, 0x08, 0x04, 0x00,
            0x00, 0x00, 0x00, 0x00, 0xff,
            0x06, 0x00, 0x42, 0x43, 0x02,
            0x00, 0x1b, 0x00, 0x03, 0x00,
            0x00, 0x00, 0x00, 0x00, 0x00,
            0x00, 0x00, 0x00]
        actual_markers = Vector{UInt8}(undef, 28)
        read!(IO, actual_markers)
        if isequal(eof_markers, actual_markers)
            @info "file have BGZF EOF Markers and is corrected!" filename = fs
        else
            @warn "file maybe truncted,please check your file!" filename = fs
        end
    end
end

function main(ARGS)
    if length(ARGS) > 1
        check_eofs.(ARGS)
    else
        check_eofs(ARGS[1])
    end
end

main(ARGS)
```

原理很简单，根据官方对SAM文件的[描述](https://samtools.github.io/hts-specs/SAMv1.pdf)，其使用了一种兼容gzip的压缩格式，但是这个格式可以通过建立索引最终达到随机访问的目的，而每个文件的末尾都存在EOF-marker, 其为
`1f 8b 08 04 00 00 00 00 00 ff 06 00 42 43 02 00 1b 00 03 00 00 00 00 00 00 00 00 00` 
的28位比特字符，我们直接读取文件末尾的28比特和该标志符对比即可检查其完整性了。

使用方法

```bash
julia check_BGZF.jl BAM_1 BAM_2 BAM_3 ...
```

同理，所有使用BGZF压缩格式的文件都可以采用这一方法检查完整，例如压缩的vcf文件`vcf.gz`。