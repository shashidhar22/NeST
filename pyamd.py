#! /usr/bin/python3
import os
import sys
import time
import shutil
import logging
import argparse
import subprocess
from pyamd.bbduk import QualCheck
from pyamd.alignment import Bwa
from pyamd.alignment import Bowtie
from pyamd.alignment import BBMap
from pyamd.alignment import Snap
from pyamd.samtools import Samtools
from pyamd.reader import Reader
from pyamd.gatk import GenAnTK
from pyamd.annotater import Annotate

def main(bbduk_path, alinger_path, smt_path, bft_path, gatk_path, 
        rone_path, rtwo_path, ref_path, adp_path, bed_path, 
        out_path, aligner):
    #Setup logging

    #Check if files are present
    if not os.path.exists(rone_path):
        raise FileNotFoundException('Forward read not found; Exiting MARs')
        sys.exit()

    if not os.path.exists(rtwo_path):
        raise FileNotFoundException('Reverse read not found; Exiting MARs')
        sys.exit()
    
    if not os.path.exists(ref_path):
        raise FileNotFoundException('Reference fasta file not found; Exiting MARs')
        sys.exit()

    if not os.path.exists(adp_path):
        raise FileNotFoundException('Adpater sequence not found; Exiting MARs')
        sys.exit()

    if not os.path.exists(out_path):
        os.mkdir(out_path)

    
    #Call Bbduk
    bbduk = QualCheck(bbduk_path, adp_path, out_path)
    rone_path, rtwo_path, bret = bbduk.bbduk(rone_path, rtwo_path)
    if bret != 0:
        raise RuntimeError('BBDuk failed to complete; Exiting MARs')

    if aligner == 'bwa':
        #Call BWA
        bwa = Bwa(alinger_path, out_path, ref_path)
        sam_path, mret = bwa.bwamem(rone_path, rtwo_path)
        if mret != 0:
            raise RuntimeError('Bwa mem failed to complete; Exiting MARs')
    elif aligner == 'bowtie2':
        #Call Bowtie2
        bowtie = Bowtie(alinger_path, out_path, ref_path)
        sam_path, mret = bowtie.bowtie(rone_path, rtwo_path)
        if mret != 0:
            raise RuntimeError('Bowtie2 failed to complete; Exiting MARs')
    elif aligner == 'snap':
        #Call Snap
        snap = Snap(alinger_path, out_path, ref_path)
        sam_path, mret = snap.snap(rone_path, rtwo_path)
        if mret != 0:
            raise RuntimeError('Snap failed to complete; Exiting MARs')
    elif aligner == 'bbmap':
        #Call Bbmap
        bbmap = BBMap(alinger_path, out_path, ref_path)
        sam_path, mret = bbmap.bbmap(rone_path, rtwo_path)
        if mret != 0:
            raise RuntimeError('BBMap failed to complete; Exitinign MARs')
    
    #Call Samtools
    varengine = Samtools(smt_path, bft_path, out_path)
    bam_path, fret = varengine.fixmate(sam_path)
    if fret != 0:
        raise RuntimeError('Samtools fixmate failed to complete; Exiting MARs')
    
    bam_path, sret = varengine.sort(sam_path)
    if sret != 0:
        raise RuntimeError('Samtools sort failed to complete; Exiting MARs')

    bcf_path, pret = varengine.pileup(ref_path, bam_path)
    if pret != 0:
        raise RuntimeError('Samtools mpileup failed to complete; Exiting MARs')

    bret = varengine.bcfindex(bcf_path)
    if bret != 0:
        raise RuntimeError('Bcftools index failed to complete; Exiting MARs')

    vcf_path, bret = varengine.bcftools(bcf_path, bed_path)
    if bret != 0:
        raise RuntimeError('Bcftools failed to complete; Exiting MARs')

    stats_path, bret = varengine.bcfstats(vcf_path, ref_path)
    if bret != 0:
        raise RuntimeError('Bcftools failed to complete; Exiting MARs')

    #Call GATK
    pic_path = 'lib/picard.jar'
    varcaller = GenAnTK(gatk_path, pic_path, out_path)
    
    add_path, aret = varcaller.picard(bam_path)
    if aret != 0:
        raise RuntimeError('GATK failed to complete; Exiting MARs')

    gvcf_path, gret = varcaller.hapCaller(add_path, ref_path)
    if gret != 0:
        raise RuntimeError('GATK failed to complete; Exiting MARs')


    #Annotate variants
    annotate = Annotate(out_path)
    annotate.iterVcf(bed_path, gvcf_path, ref_path, 'gatk')
    annotate.iterVcf(bed_path, vcf_path, ref_path, 'samtools')
   
if __name__ == '__main__':

    bbduk_def = shutil.which("bbduk.sh")
    bwa_def = shutil.which("bwa")
    bowtie_def = shutil.which("bowtie2")
    snap_def = shutil.which("snap-aligner")
    smt_def = shutil.which("samtools")
    bft_def = shutil.which("bcftools")
    gatk_def = shutil.which("GenomeAnalysisTK.jar")

    #Get arguments
    parser = argparse.ArgumentParser(prog='kookaburra')
    parser.add_argument('-1', '--fwd', dest='rone_path', type=str, 
                        help='Path to forward reads fastq', required=True)
    parser.add_argument('-2', '--rev', dest='rtwo_path', type=str, 
                        help='Path to reverse reads fastq', required=True)
    parser.add_argument('-r', '--ref', dest='ref_path', type=str, 
                        help='Path to Reference fasta file', required=True)
    parser.add_argument('-a', '--adapter', dest='adp_path', type=str, 
                        help='Path to Adpater fasta file', required=True)
    parser.add_argument('-b', '--bed', dest='bed_path', type=str, 
                        help='Path to Bed file for MDR regions', required=True)
    parser.add_argument('-o', '--outpath', dest='out_path', type=str, 
                        help='Path where all outputs will be stored', required=True)
    parser.add_argument('-m', '--mapper', dest='aligner', type=str,
                        choices=['bowtie2', 'bwa', 'bbmap', 'snap'], 
                        default='bwa', help='The aligner to used by MARs')
    parser.add_argument('--bbduk', dest='bbduk_path', type=str, default=bbduk_def,
                        help='Path to BBduk executable')
    parser.add_argument('--aligner', dest='aligner_path', type=str, default=bwa_def,
                        help='Path to aligner executable')
    parser.add_argument('--samtools', dest='smt_path', type=str, default=smt_def,
                        help='Path to Samtools executable')
    parser.add_argument('--gatk', dest='gatk_path', type=str, default=gatk_def,
                        help='Path to GATK executable')
    parser.add_argument('--bcftools', dest='bft_path', type=str, default=bft_def,
                        help='Path to Bcftools executable')

    

    
    args = parser.parse_args()
    main(args.bbduk_path, args.aligner_path, args.smt_path, args.bft_path, args.gatk_path, 
        args.rone_path, args.rtwo_path, args.ref_path, args.adp_path, args.bed_path, 
        args.out_path, args.aligner)
