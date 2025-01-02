#!/bin/bash

ROI="../../data/opsin_region.bed"
WINDOW=1000000
REF="Homo_sapiens.GRCh38.dna.chromosome.X.fa"
URL="https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.X.fa.gz"
OUT=grch38.chrX.roi.fa

wget -O "${REF}.gz" ${URL}
gunzip "${REF}.gz"
sed -i -e 's/^>/>chr/' "${REF}"
samtools faidx ${REF}
cat ${REF}.fai | awk 'BEGIN{FS="\t"; OFS=FS} {print $1,"0",$2}' > ${REF}.bed
bedtools slop -i ${ROI} -g ${REF}.fai -b $WINDOW > roi_padded.bed
bedtools subtract -a ${REF}.bed -b roi_padded.bed > mask_region.bed
bedtools maskfasta -fi ${REF} -bed mask_region.bed -fo ${OUT}
bgzip ${OUT}
samtools faidx ${OUT}.gz

rm ${REF}
rm ${REF}.bed
rm ${REF}.fai
rm mask_region.bed
rm roi_padded.bed