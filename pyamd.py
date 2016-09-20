#! /usr/bin/python3
import os
import sys
import time
import shutil
import logging
import argparse
import subprocess
from pyamd.bbduk import QualCheck
from pyamd.bwa import Bwa
from pyamd.samtools import Samtools
from pyamd.reader import Reader
from pyamd.gatk import GenAnTK

def main(bbduk_path, bwa_path, smt_path, bft_path, gatk_path, 
        rone_path, rtwo_path, ref_path, adp_path, bed_path, out_path):
    #Setup logging

    #Check if files are present
    if not os.path.exists(rone_path):
        raise FileNotFoundException('Forward read not found; Exiting pyAMD')
        sys.exit()

    if not os.path.exists(rtwo_path):
        raise FileNotFoundException('Reverse read not found; Exiting pyAMD')
        sys.exit()
    
    if not os.path.exists(ref_path):
        raise FileNotFoundException('Reference fasta file not found; Exiting pyAMD')
        sys.exit()

    if not os.path.exists(adp_path):
        raise FileNotFoundException('Adpater sequence not found; Exiting pyAMD')
        sys.exit()

    if not os.path.exists(out_path):
        os.mkdir(out_path)

    
    #Call Bbduk
    bbduk = QualCheck(bbduk_path, adp_path, out_path)
    rone_path, rtwo_path, bret = bbduk.bbduk(rone_path, rtwo_path)
    if bret != 0:
        raise RuntimeError('BBDuk failed to complete; Exiting pyAMD')
    #Call BWA
    bwa = Bwa(bwa_path, out_path, ref_path)
    sam_path, mret = bwa.bwamem(rone_path, rtwo_path)
    if mret != 0:
        raise RuntimeError('Bwa mem failed to complete; Exiting pyAMD')

    #Call Samtools
    varengine = Samtools(smt_path, bft_path, out_path)
    bam_path, fret = varengine.fixmate(sam_path)
    if fret != 0:
        raise RuntimeError('Samtools fixmate failed to complete; Exiting pyAMD')
    
    bam_path, sret = varengine.sort(sam_path)
    if sret != 0:
        raise RuntimeError('Samtools sort failed to complete; Exiting pyAMD')

    bcf_path, pret = varengine.pileup(ref_path, bam_path)
    if pret != 0:
        raise RuntimeError('Samtools mpileup failed to complete; Exiting pyAMD')

    bret = varengine.bcfindex(bcf_path)
    if bret != 0:
        raise RuntimeError('Bcftools index failed to complete; Exiting pyAMD')

    vcf_path, bret = varengine.bcftools(bcf_path, bed_path)
    if bret != 0:
        raise RuntimeError('Bcftools failed to complete; Exiting pyAMD')

    stats_path, bret = varengine.bcfstats(vcf_path, ref_path)
    if bret != 0:
        raise RuntimeError('Bcftools failed to complete; Exiting pyAMD')

    #Call GATK
    pic_path = 'lib/picard.jar'
    varcaller = GenAnTK(gatk_path, pic_path, out_path)
    
    add_path, aret = varcaller.picard(bam_path)
    if aret != 0:
        raise RuntimeError('GATK failed to complete; Exiting pyAMD')

    gvcf_path, gret = varcaller.hapCaller(add_path, ref_path)
    if gret != 0:
        raise RuntimeError('GATK failed to complete; Exiting pyAMD')
   
    #Extract variants in bed region 
    reader = Reader()
    ovars = reader.extractVars(vcf_path, ref_path)
    obed_path = '{0}/variants.bed'.format(out_path)
    obed = open(obed_path, 'w')
    obed.write('Chrom\tPos\tRef\tAlt\tCodonPos\trefCodon\taltCodon\n')
    for var in ovars:
        obed.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(var.chrom, var.pos, var.ref, var.alt, var.codonPos, var.refCodon, var.altCodon))
    obed.close()

    ovars = reader.extractVars(gvcf_path, ref_path)
    obed_path = '{0}/variants_gatk.bed'.format(out_path)
    obed = open(obed_path, 'w')
    obed.write('Chrom\tPos\tRef\tAlt\tCodonPos\trefCodon\taltCodon\n')
    for var in ovars:
        obed.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(var.chrom, var.pos, var.ref, var.alt, var.codonPos, var.refCodon, var.altCodon))
    obed.close()

if __name__ == '__main__':
    #Set default tool paths
#def main(bbduk_path, bwa_path, smt_path, bft_path, rone_path, 
#        rtwo_path, ref_path, adp_path, out_path):

    bbduk_def = shutil.which("bbduk.sh")
    bwa_def = shutil.which("bwa")
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
    parser.add_argument('--bbduk', dest='bbduk_path', type=str, default=bbduk_def,
                        help='Path to BBduk executable')
    parser.add_argument('--bwa', dest='bwa_path', type=str, default=bwa_def,
                        help='Path to Bwa executable')
    parser.add_argument('--samtools', dest='smt_path', type=str, default=smt_def,
                        help='Path to Samtools executable')
    parser.add_argument('--gatk', dest='gatk_path', type=str, default=gatk_def,
                        help='Path to GATK executable')
    parser.add_argument('--bcftools', dest='bft_path', type=str, default=bft_def,
                        help='Path to Bcftools executable')

    

    
    args = parser.parse_args()
    main(args.bbduk_path, args.bwa_path, args.smt_path, args.bft_path, args.gatk_path,
        args.rone_path, args.rtwo_path, args.ref_path, args.adp_path, args.bed_path, args.out_path)
