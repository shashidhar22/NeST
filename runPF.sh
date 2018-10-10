#!/bin/bash
mkdir -p local
time python3 nest.py -i fq/MaRS_test/SRR_Acc_List.txt -a ref/pfalciparum/adapters.fa -r ref/pfalciparum/mdr.fa -o local/MaRS_test -b ref/pfalciparum/mdr.bed -m bowtie2 --varofint ref/pfalciparum/Reportable_SNPs.csv
