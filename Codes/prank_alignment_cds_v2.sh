#!/bin/sh

##this takes a list of file names without the extension
##the file name can be a protein name, ortholog name identifier or such (written by a previous program)
##it aligns all coding sequences within the .fas file identified by the name and writes to a .aln file with same name
##plus, it logs its actions

fn_oi=$1
wt_dir=$2
wd_oi=$3

echo $fn_oi

cds_dir="${wt_dir}${wd_oi}/${wd_oi}_cds_alignments/"

echo "PRANK was called as prank -d=${cds_dir}${fn_oi}_cds.fas -o=${cds_dir}${fn_oi} -translate -F" >> ${cds_dir}prank_translate_alignment.log

prank -d=${cds_dir}${fn_oi}_cds.fas -o=${cds_dir}${fn_oi} -translate -F >> ${cds_dir}prank_translate_alignment.log
