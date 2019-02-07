#!/bin/bash
python3 prepSubmission.py
mkdir -p local
time python3 nest.py -i fq/ColmanEtAl/ -a ref/mtuberculosis/adapters.fa -r ref/mtuberculosis/mtuberculosis.fa -o local/ColmanEtAl/ -b ref/mtuberculosis/tbmdr.bed -m bowtie2 --varofint ref/mtuberculosis/coleman.csv
