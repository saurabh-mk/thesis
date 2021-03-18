#!/bin/sh

## this takes the name of an multiple alignment file and cleans it using gblocks

fn_oi=$1
wt_dir=$2
wd_oi=$3

cds_dir="${wd_oi}_cds_alignments/"

cd ${wt_dir}${wd_oi}/

echo "Gblocks was called as Gblocks ${cds_dir}${fn_oi}.best.nuc.fas -t=c -b1=12 -b2=12 -b3=2 -b4=10 -b5=n -e=aln" >> ${cds_dir}gblocks_alignment.log

Gblocks ${cds_dir}${fn_oi}.best.nuc.fas -t=c -b1=12 -b2=12 -b3=2 -b4=10 -b5=n -e=.aln >> ${cds_dir}gblocks_alignment.log
