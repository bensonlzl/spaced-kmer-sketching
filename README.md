# Spaced K-mer Sketching

This codebase provides a (somewhat) performant set of functions for K-mer sketching with spaced seeds in C++.

## What are K-mers?

K-mers are contiguous substrings of a fixed length, usually used in the context of genomic or proteomic data.
For example, the nucleotide string `ACCGTAAATTCGA` consists of the 5-mers

- `ACCGT`
- `CCGTA`
- `CGTAA`
- `GTAAA`
- `TAAAT`
- `AAATT`
- `AATTC`
- `ATTCG`
- `TTCGA`

## What are spaced seeds?

Spaced seeds take k-mers and insert wildcard positions in the k-mer,
positions where we simply ignore the nucleotide.

For example, if we are looking at a standard 8-mer `ACGTACGT` 
from `AAACGTACGTTT`, we can view it as choosing a window length of 8
and using the seed `11111111`

```
AAACGTACGTTT
  ||||||||  
..ACGTACGT..
```

Now if we use a different seed, say `11001011`, we get the 5-mer `ACAGT` 

```
AAACGTACGTTT
  ||  | ||  
..AC..A.GT..
```

## What is Sketching?

