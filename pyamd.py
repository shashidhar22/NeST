#! /usr/bin/python3
import os
import sys
import glob
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
from pyamd.gatk import Picard
from pyamd.annotater2 import Annotate
from pyamd.kestrel import kes_runner
from pyamd.filter import filterer
from pyamd.summarize import Summary
from pyamd.prepinputs import Prepper

def main(bbduk_path, alinger_path, smt_path, bft_path, gatk_path,
        rone_path, rtwo_path, ref_path, adp_path, bed_path,
        out_path, aligner,kes_path, kan_path, pic_path, sam_name):
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

    #Fix mate information, sort files and add read groups
    varengine = Samtools(smt_path, bft_path, out_path)
    bam_path, fret = varengine.fixmate(sam_path)
    if fret != 0:
        raise RuntimeError('Samtools fixmate failed to complete; Exiting MARs')

    bam_path, sret = varengine.sort(sam_path)
    if sret != 0:
        raise RuntimeError('Samtools sort failed to complete; Exiting MARs')

    rgadder = Picard(pic_path, out_path)
    bam_path, aret = rgadder.picard(bam_path, sam_name)
    if aret != 0:
        raise RuntimeError('GATK failed to complete; Exiting MARs')

    bcf_path, pret = varengine.pileup(ref_path, bam_path)
    if pret != 0:
        raise RuntimeError('Samtools mpileup failed to complete; Exiting MARs')

    bret = varengine.bcfindex(bcf_path)
    if bret != 0:
        raise RuntimeError('Bcftools index failed to complete; Exiting MARs')

    vcf_path, bret = varengine.bcftools(bcf_path, bed_path, sam_name)
    if bret != 0:
        raise RuntimeError('Bcftools failed to complete; Exiting MARs')

    stats_path, bret = varengine.bcfstats(vcf_path, ref_path)
    if bret != 0:
        raise RuntimeError('Bcftools failed to complete; Exiting MARs')

    #Call GATK
    #pic_path = 'lib/picard.jar'
    varcaller = GenAnTK(gatk_path, out_path)


    gvcf_path, gret = varcaller.hapCaller(bam_path, ref_path, sam_name)
    if gret != 0:
        raise RuntimeError('GATK failed to complete; Exiting MARs')



    merged_vcf = filterer(gvcf_path, vcf_path, sam_name, out_path)
    annotate = Annotate(out_path)
    annotate.iterVcf(bed_path, merged_vcf, sam_name, ref_path, 'merged')
    annotate.iterVcf(bed_path, gvcf_path, sam_name, ref_path, 'gatk')
    annotate.iterVcf(bed_path, vcf_path , sam_name, ref_path, 'samtools')
    return

def marsBatch(bbduk_path, aligner_path, smt_path, bft_path, gatk_path,
              inp_path, ref_path, adp_path, bed_path,
              out_dir, aligner, kes_path, kan_path, pic_path, voi_path):
    #Setup loggers
    logger = logging.getLogger('MaRS Batch')
    logger.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    formatter = logging.Formatter('%(asctime)s : %(name)s : %(levelname)s : %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    #Create output directories
    
    if not os.path.exists(os.path.abspath(out_dir)):
        os.mkdir(os.path.abspath(out_dir))
    logger.info('Output will be stored in {0}'.format(os.path.abspath(out_dir)))        
    
    #Prepare input configuration
    prep = Prepper(inp_path)
    config = prep.prepInputs()
    for samples in config:
        sample_name = config[samples].sample
        rone_path = config[samples].files[0]
        rtwo_path = config[samples].files[1]
        out_path = '{0}/{1}'.format(os.path.abspath(out_dir), sample_name)
        if not os.path.exists(out_path):
            os.mkdir(out_path)
        main(bbduk_path, aligner_path, smt_path, bft_path, gatk_path,
             rone_path, rtwo_path, ref_path, adp_path, bed_path,
             out_path, aligner, kes_path, kan_path, pic_path, sample_name)
    summary = Summary(ref_path, bed_path, voi_path, out_dir)
    summary.getVarOfInt(out_dir)


if __name__ == '__main__':

    def_path = "{0}/lib".format(os.path.abspath(os.path.dirname(os.path.realpath(__file__))))
    bbduk_def = "{0}/bbmap/bbduk.sh".format(def_path)
    bwa_def = "{0}/bwa-0.7.12/bwa".format(def_path)
    bowtie_def = "{0}/bowtie2-2.3.0/bowtie2".format(def_path)
    snap_def = "{0}/snap/snap-aligner".format(def_path)
    smt_def = "{0}/samtools-1.3.1/samtools".format(def_path)
    bft_def = "{0}/bcftools-1.3.1/bcftools".format(def_path)
    gatk_def = "{0}/GenomeAnalysisTK.jar".format(def_path)
    pic_def = "{0}/picard.jar".format(def_path)
    aligner_def = {'bwa' : bwa_def, 'snap' : snap_def, 'bowtie2': bowtie_def}
    #Get arguments
    parser = argparse.ArgumentParser(prog='kookaburra')
    parser.add_argument('-i', '--inp_path', type=str,
                        help='Path to input directory (Specify only for batch mode)')
    parser.add_argument('-1', '--fwd', dest='rone_path', type=str,
                        help='Path to forward reads fastq', )
    parser.add_argument('-2', '--rev', dest='rtwo_path', type=str,
                        help='Path to reverse reads fastq')
    parser.add_argument('-r', '--ref', dest='ref_path', type=str,
                        help='Path to Reference fasta file', required=True)
    parser.add_argument('-a', '--adapter', dest='adp_path', type=str,
                        help='Path to Adpater fasta file', required=True)
    parser.add_argument('-b', '--bed', dest='bed_path', type=str,
                        help='Path to Bed file for MDR regions', required=True)
    parser.add_argument('-o', '--outpath', dest='out_path', type=str,
                        help='Path where all outputs will be stored', required=True)
    parser.add_argument('-n', '--sam_name', dest='sam_name', type=str,
                        help='Sample name', default=None)
    parser.add_argument('-m', '--mapper', dest='aligner', type=str,
                        choices=['bowtie2', 'bwa', 'bbmap', 'snap'],
                        default='bwa', help='The aligner to used by MARs')
    parser.add_argument('--bbduk', dest='bbduk_path', type=str, default=bbduk_def,
                        help='Path to BBduk executable')
    parser.add_argument('--aligner', dest='aligner_path', type=str, default=None,
                        help='Path to aligner executable')
    parser.add_argument('--samtools', dest='smt_path', type=str, default=smt_def,
                        help='Path to Samtools executable')
    parser.add_argument('--gatk', dest='gatk_path', type=str, default=gatk_def,
                        help='Path to GATK executable')
    parser.add_argument('--bcftools', dest='bft_path', type=str, default=bft_def,
                        help='Path to Bcftools executable')
    parser.add_argument('--picard', dest='pic_path', type=str, default=pic_def,
                        help='Path to Bcftools executable')
    parser.add_argument('--kestrel', dest='kes_path', type=str, default=bft_def,
                        help='Path to Kestrel executable')
    parser.add_argument('--kanalyze', dest='kan_path', type=str, default=bft_def,
                        help='Path to Kanalyze executable')

    parser.add_argument('--varofint', dest='voi_path', type=str, default=bft_def,
                        help='Path to variant of interest')


    
    args = parser.parse_args()
    if args.aligner_path == None:
        args.aligner_path = aligner_def[args.aligner]
 
    if not os.path.exists(args.out_path):
        os.mkdir(args.out_path)

    if args.inp_path == None and args.rone_path != None:
        if args.sam_name == None:
            sam_name = os.path.splitext(os.path.basename(args.rone_path))[0]
        else:
            sam_name = args.sam_name
        main(args.bbduk_path, args.aligner_path, args.smt_path, args.bft_path, args.gatk_path,
            args.rone_path, args.rtwo_path, args.ref_path, args.adp_path, args.bed_path,
            args.out_path, args.aligner, args.kes_path, args.kan_path, args.pic_path, sam_name)
    elif args.inp_path != None and args.rone_path == None:
        marsBatch(args.bbduk_path, args.aligner_path, args.smt_path, args.bft_path, args.gatk_path,
            args.inp_path, args.ref_path, args.adp_path, args.bed_path,
            args.out_path, args.aligner, args.kes_path, args.kan_path, args.pic_path, args.voi_path)
