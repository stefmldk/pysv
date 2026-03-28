# PySV
A fast Python-based, long-read structural variant caller for Linux

## Installation
Simply download pysv.py and pysv_env.yml. Use pysv_env.yml to create a Conda environment with the required dependencies: Assuming mamba is installed, run the command
```
mamba env create -f pysv_env.yml
```

## Quick start
To get a quick overview of command parameters run
```
python pysv.py -h
```
To run PySV using four cores for BAM file searching, use the command:
```
python pysv.py -b <bam_file> -c4
```
To run PySV recursively on all bam files in a folder, use the command:
```
python pysv.py -d <folder> -c4
```
## Overview
Pysv is designed to call large structural variants (SVs) by analyzing chimeric reads. Chimeric reads are reads that contain multiple supplementary alignments that map to different places in the reference genome. It does not in the current version analyse CIGAR sequences in order to call SVs that are small enough to be contained within individual alignments. Other tools may be used for that purpose.
