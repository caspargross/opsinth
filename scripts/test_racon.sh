#!/bin/bash


./lib/racon \
    -t 1 \
    tests/out/crispr/opsinth.reads.fastq \
    tests/out/crispr/opsinth.denovo.unsorted.sam \
    tests/out/crispr/opsinth.denovo.fasta \
    > tests/out/crispr/opsinth.denovo.corrected.fasta

