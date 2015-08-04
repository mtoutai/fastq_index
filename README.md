THIS SOFTWARE IS DISTRIBUTED UNDER
CC0 1.0 Universal
Public Domain Dedication
http://creativecommons.org/publicdomain/zero/1.0/

Developed by:
TOTAI MITSUYAMA, PhD.
Biotechnology Research Institute for Drug Discovery
National Institute of Advanced Industrial Science and Technology (AIST)



INTRODUCTION
In many cases of Next-Generation Sequencer (NGS) data analyses, FASTQ files are split into several parts in order to parallelize jobs and gain the best performance out of a cluster computing environment. However, splitting large FASTQ files consumes the precious disk space. Sometimes, such a situation is not acceptable. FASTQ_INDEX provides a simple alternative to splitting files. It creates an index file for a given FASTQ file which enables random read access to a subset of sequences in a FASTQ file similar to accessing a split file of a FASTQ file.
A typical use case of aligning short reads with FASTQ_INDEX looks like this:

    > fastq_index test.fastq
    > ls -1
    test.fastq
    test.fastq.index
    > fastq_index_cat n 3 test.fastq | lastal Q1 hg38.fa - > result.maf


This software is a set of following tools:
- fastq_index         creates an index file for a given FASTQ file
- fastq_index_cat     writes specified block of a FASTQ file to STDOUT
- fastq_index_dump    writes contents of an index file to STDOUT (debug purpose)
- fastq_index_size    writes the number of blocks of an index file to STDOUT


INSTALLATION

A typical installation procedure is to simply issue a make command.

    > make

Then, you can manually copy executable to where you want to store them.


