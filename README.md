# PYSV
A fast Python-based, long-read structural variant caller for Linux

## Installation
Simply download pysv.py and

Pysv is designed to call relatively large structural variants (SVs), like large del/dups, inversions or translocations. It does not analyse CIGAR sequences and therefore will not call SVs small enough to be contained in these. Other tools exist that can call SVs contained in CIGAR sequences.

In the first part of a run, the given BAM file is searched for reads containing supplementary alignments. This part can be accelerated via multiprocessing.

Next, pairs of consecutive alignments (consecutive on the reads, not on the reference) are created and clusters of pairs that share breakpoints are created. A maximum breakpoint wobble may be defined for this process. Default is +/- 20 bases.

The resulting clusters and their associated brekpoints may at this point be used for manual examination in a suitable genome browser like IGV. However, PYSV 

bases its algorithm on reads containing supplementary alignments. It does not in this version examine the CIGAR sequences and therefore it is primarily relatively large structural variants (SVs) that pysv detects
