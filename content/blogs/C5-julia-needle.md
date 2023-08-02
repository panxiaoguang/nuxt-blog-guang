---
title: Julia语言编写Needleman Wunsch全局比对算法
date: 23th Mar 2023
description: Julia语言编写Needleman Wunsch全局比对算法.
image: /blogs-img/julia.webp
alt: Julia语言编写Needleman Wunsch全局比对算法
ogImage: /blog-img/julia.webp
tags: ['julia']
published: true
---

输入是两个字符串，输出是对齐后的两个字符串。
```r
function needleman_wunsch(seq1::String, seq2::String; gap_penalty::Int64=-1, match_score::Int64=1, mismatch_penalty::Int64=-1)
    # Check if the sequences are valid
    if !all(isuppercase(i) for i in seq1) || !all(isuppercase(j) for j in seq2)
        throw(ArgumentError("The sequences must contain only upper-case letters."))
    end
    # Calculate the length of the two sequences
    len1, len2 = length(seq1), length(seq2)
    # Initialize the matrix
    matrix = zeros(Int64, len1 + 1, len2 + 1)
    # Initialize the first column
    first_col = @views matrix[2:end, 1]
    first_col[:] .= [i * gap_penalty for i in 1:(length(seq1))]
    # Initialize the first row
    first_row = @views matrix[1, 2:end]
    first_row[:] .= [j*gap_penalty for j in 1:(length(seq2))]
    # Fill the matrix
    for i in 2:(len1+1)
        for j in 2:(len2+1)
            match = seq1[i-1] == seq2[j-1] ? matrix[CartesianIndex(i - 1, j - 1)] + match_score : matrix[CartesianIndex(i - 1, j - 1)] + mismatch_penalty
            delete = matrix[CartesianIndex(i - 1, j)] + gap_penalty
            insert = matrix[CartesianIndex(i, j - 1)] + gap_penalty
            matrix[CartesianIndex(i, j)] = max(match, delete, insert)
        end
    end
    # Initialize the aligned sequences
    aligned_seq1 = ""
    aligned_seq2 = ""
    # Initialize the current index
    i = len1 + 1
    j = len2 + 1
    # Traceback
    while i > 1 || j > 1
        if i > 1 && j > 1 && matrix[i, j] == matrix[i-1, j-1] + (seq1[i-1] == seq2[j-1] ? match_score : mismatch_penalty)
            aligned_seq1 = string(seq1[i-1], aligned_seq1)
            aligned_seq2 = string(seq2[j-1], aligned_seq2)
            i -= 1
            j -= 1
        elseif i > 1 && matrix[i, j] == matrix[i-1, j] + gap_penalty
            aligned_seq1 = string(seq1[i-1], aligned_seq1)
            aligned_seq2 = string('-', aligned_seq2)
            i -= 1
        else
            aligned_seq1 = string('-', aligned_seq1)
            aligned_seq2 = string(seq2[j-1], aligned_seq2)
            j -= 1
        end
    end
    (aligned_seq1, aligned_seq2)
end
```

例如：

```r
@time needleman_wunsch("AGTACGCAGTCAGGTGACGTA","AGTACCCAGTCAGGGACGTA")
  0.000017 seconds (129 allocations: 8.109 KiB)
("AGTACGCAGTCAGGTGACGTA", "AGTACCCAGTCAGG-GACGTA")
```
整个算法主要是矩阵操作。