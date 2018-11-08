#! /usr/bin/python3
import os
import sys
import glob
import time
import shutil
import logging
import argparse
import subprocess
import pandas as pd
from pathlib import Path
from itertools import repeat
from multiprocessing import Pool
from nest.bbduk import QualCheck
from nest.alignment import Bwa
from nest.alignment import Bowtie
from nest.alignment import BBMap
from nest.alignment import Snap
from nest.samtools import Samtools
from nest.gatk import GenAnTK
from nest.gatk import Picard
from nest.kestrel import KestrelVar
#from nest.annotater import Annotate
from nest.kestrel import kes_runner
from nest.summarize import Summary
from nest.prepinputs import Prepper
from nest.parsers.vcf import Vcf

def main(arguments):
    bbduk_path = arguments[0]
    alinger_path = arguments[1]
    smt_path = arguments[2]
    bft_path = arguments[3]
    gatk_path = arguments[4]
    rone_path = arguments[5]
    rtwo_path = arguments[6]
    ref_path = arguments[7]
    adp_path = arguments[8]
    bed_path = arguments[9]
    out_dir = arguments[10]
    aligner = arguments[11]
    pic_path = arguments[12]
    sam_name = arguments[13]
    voi_path = arguments[14]
    java_path = arguments[15]
    #Setup logging
    #Get logger for main method
    main_logger = logging.getLogger('NeST.{0}'.format(sam_name))

    #Check if files are present
    out_path = '{0}/{1}'.format(os.path.abspath(out_dir), sam_name)
    if not os.path.exists(out_path):
        os.mkdir(out_path)


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

    #Create completion folder
    completion_path = '{0}/completion'.format(out_path)
    if not os.path.exists(completion_path):
        os.mkdir(completion_path)

    #Call Bbduk
    main_logger.debug('Running BBDuk')
    if os.path.exists('{0}/bbduk.rt'.format(completion_path)):
        brone = os.path.splitext(os.path.basename(rone_path))[0]
        brtwo = os.path.splitext(os.path.basename(rtwo_path))[0]
        rone_path = '{0}/{1}/{2}_cleaned.fq'.format(out_path, 'CleanedFastq', brone)
        rtwo_path = '{0}/{1}/{2}_cleaned.fq'.format(out_path, 'CleanedFastq', brtwo)
        main_logger.debug('Skipping BBDuk')
        bret = 0
    else:
        bbduk = QualCheck(bbduk_path, adp_path, out_path, java_path)
        rone_path, rtwo_path, bret = bbduk.bbduk(rone_path, rtwo_path)
        if bret == 0:
            Path('{0}/bbduk.rt'.format(completion_path)).touch()
    if bret != 0:
        raise RuntimeError('BBDuk failed to complete; Exiting MARs')
    else:
        main_logger.debug('BBDuk completed')

    if aligner == 'bwa':
        #Call BWA
        main_logger.debug('Running BWA')
        if os.path.exists('{0}/align.rt'.format(completion_path)):
            sam_path = '{0}/alignments/output.sam'.format(out_path)
            mret = 0
            main_logger.debug('Skipping BWA')
        else:
            bwa = Bwa(alinger_path, out_path, ref_path)
            sam_path, mret = bwa.bwamem(rone_path, rtwo_path)
            if mret == 0:
                Path('{0}/align.rt'.format(completion_path)).touch()
        if mret != 0:
            raise RuntimeError('Bwa mem failed to complete; Exiting MARs')
        else:
            main_logger.debug('BWA completed')

    elif aligner == 'bowtie2':
        #Call Bowtie2
        main_logger.debug('Running Bowtie2')
        if os.path.exists('{0}/aling.rt'.format(completion_path)):
            sam_path = '{0}/alignments/output.sam'.format(out_path)
            mret = 0
            main_logger.debug('Skipping Bowtie2')
        else:
            bowtie = Bowtie(alinger_path, out_path, ref_path)
            sam_path, mret = bowtie.bowtie(rone_path, rtwo_path)
            if mret == 0:
                Path('{0}/align.rt'.format(completion_path)).touch()
        if mret != 0:
            raise RuntimeError('Bowtie2 failed to complete; Exiting MARs')
        else:
            main_logger.debug('Bowtie2 completed')

    elif aligner == 'snap':
        #Call Snap
        main_logger.debug('Running Snap')
        snap = Snap(alinger_path, out_path, ref_path)
        sam_path, mret = snap.snap(rone_path, rtwo_path)
        if mret != 0:
            raise RuntimeError('Snap failed to complete; Exiting MARs')
        else:
            main_logger.debug('Snap completed')

    elif aligner == 'bbmap':
        #Call Bbmap
        main_logger.debug('Running BBMap')
        if os.path.exists('{0}/aling.rt'.format(completion_path)):
            sam_path = '{0}/alignments/output.sam'.format(out_path)
            mret = 0
        else:
            bbmap = BBMap(alinger_path, out_path, ref_path)
            sam_path, mret = bbmap.bbmap(rone_path, rtwo_path)
            if mret == 0:
                Path('{0}/align.rt'.format(completion_path)).touch()
        if mret != 0:
            raise RuntimeError('BBMap failed to complete; Exitinign MARs')
        else:
            main_logger.debug('BBMap completed')


    #Fix mate information, sort files and add read groups
    varengine = Samtools(smt_path, bft_path, out_path)
    if os.path.exists('{0}/fixmate.rt'.format(completion_path)):
        base = os.path.splitext(os.path.basename(sam_path))[0]
        bam_path = '{0}/alignments/{1}_FM.bam'.format(out_path, base)
        fret = 0
        main_logger.debug('Skipping fixmate')
    else:
        bam_path, fret = varengine.fixmate(sam_path)
        if fret == 0:
            Path('{0}/fixmate.rt'.format(completion_path)).touch()
    main_logger.debug('Running Samtools fixmate')
    if fret != 0:
        raise RuntimeError('Samtools fixmate failed to complete; Exiting MARs')
    else:
        main_logger.debug('Samtools fixmate completed')

    if os.path.exists('{0}/sort.rt'.format(completion_path)):
        base = os.path.splitext(os.path.basename(bam_path))[0]
        bam_path = '{0}/alignments/{1}_SR.bam'.format(out_path, base)
        sret = 0
        main_logger.debug('Skipping sort')
    else:
        bam_path, sret = varengine.sort(bam_path)
        if sret == 0:
            Path('{0}/sort.rt'.format(completion_path)).touch()
    main_logger.debug('Running Samtools sort')
    if sret != 0:
        raise RuntimeError('Samtools sort failed to complete; Exiting MARs')
    else:
        main_logger.debug('Samtools sort completed')

    if os.path.exists('{0}/dedup.rt'.format(completion_path)):
        base = os.path.splitext(os.path.basename(bam_path))[0]
        bam_path = '{0}/alignments/{1}_DD.bam'.format(out_path, base)
        dret = 0
        main_logger.debug('Skipping Dedup')
    else:
        bam_path, dret = varengine.dedup(bam_path)
        if dret == 0:
            Path('{0}/dedup.rt'.format(completion_path)).touch()
    main_logger.debug('Running Samtools dedup')
    if sret != 0:
        raise RuntimeError('Samtools dedup failed to complete; Exiting MARs')
    else:
        main_logger.debug('Samtools dedup completed')


    rgadder = Picard(java_path, pic_path, out_path)
    if os.path.exists('{0}/readgroup.rt'.format(completion_path)):
        base = os.path.splitext(os.path.basename(bam_path))[0]
        bam_path = '{0}/alignments/{1}_RG.bam'.format(out_path, base)
        aret = 0
        main_logger.debug('Skipping add read group')
    else:
        bam_path, aret = rgadder.picard(bam_path, sam_name)
        main_logger.debug('Running Picard AddOrReplaceReadGroups')
        if aret == 0:
            Path('{0}/readgroup.rt'.format(completion_path)).touch()
    if aret != 0:
        raise RuntimeError('Picard AddOrReplaceReadGroups failed to complete; Exiting MARs')
    else:
        main_logger.debug('Picard AddOrReplaceReadGroups completed')

    #Run samtools mpileup, bcftools index, call and stats to generate VCF files
    if os.path.exists('{0}/pileup.rt'.format(completion_path)):
        bcf_path = '{0}/{1}_variants.bcf'.format(out_path, sam_name)
        pret = 0
        main_logger.debug('Skipping Pileup')
    else:
        bcf_path, pret = varengine.pileup(ref_path, bam_path, sam_name)
        main_logger.debug('Running Samtools mpileup')
        if pret == 0:
            Path('{0}/pileup.rt'.format(completion_path)).touch()
    if pret != 0:
        raise RuntimeError('Samtools mpileup failed to complete; Exiting MARs')
    else:
        main_logger.debug('Samtools mpileup completed')

    if os.path.exists('{0}/bcfindex.rt'.format(completion_path)):
        bret = 0
        main_logger.debug('Skipping Bcfindex')
    else:
        bret = varengine.bcfindex(bcf_path)
        main_logger.debug('Running Bcftools index')
        if bret ==0 :
            Path('{0}/bcfindex.rt'.format(completion_path)).touch()
    if bret != 0:
        raise RuntimeError('Bcftools index failed to complete; Exiting MARs')
    else:
        main_logger.debug('Bcftools index completed')

    if os.path.exists('{0}/bcfcall.rt'.format(completion_path)):
        vcf_path = '{0}/{1}_variants_samtools.vcf'.format(out_path, sam_name)
        bret = 0
        main_logger.debug('Skipping bcfcall')
    else:
        vcf_path, bret = varengine.bcftools(bcf_path, sam_name)
        main_logger.debug('Running Bcftools call')
        if bret == 0:
            Path('{0}/bcfcall.rt'.format(completion_path)).touch()

    if bret != 0:
        raise RuntimeError('Bcftools call failed to complete; Exiting MARs')
    else:
        main_logger.debug('Bcftools call completed')


    #Call GATK HaplotypeCaller to generate VCF files
    varcaller = GenAnTK(gatk_path, out_path, java_path, pic_path)
    main_logger.debug('Running GATK HaplotypeCaller')
    if os.path.exists('{0}/gatk.rt'.format(completion_path)):
        gvcf_path = '{0}/{1}_variants_gatk.vcf'.format(out_path, sam_name)
        gret = 0
        main_logger.debug('Skipping GATK')
    else:
        gvcf_path, gret = varcaller.hapCaller(bam_path, ref_path, sam_name)
        if gret == 0:
            Path('{0}/gatk.rt'.format(completion_path)).touch()
    if gret != 0:
        raise RuntimeError('GATK HaplotypeCaller failed to complete; Exiting MARs')
    else:
        main_logger.debug('GATK HaplotypeCaller stats completed')

    #Call Kestrel to generate VCF files
    #kestrel_path = 'lib/kestrel/kestrel.jar'
    #kanalyze_path = 'lib/kestrel/kanalyze.jar'
    #varcaller = KestrelVar(rone_path, rtwo_path, ref_path, kanalyze_path,
    #                        kestrel_path, out_path)
    #varcaller = GenAnTK(gatk_path, out_path, java_path)
    #main_logger.debug('Running Kestrel')
    #if os.path.exists('{0}/kestrel.rt'.format(completion_path)):
    #    kvcf_path = '{0}/vairants_kes.vcf'.format(out_path)
    #    kret = 0
    #    main_logger.debug('Skipping Kestrel')
    #else:
    #    kvcf_path, kret = varcaller.run_kestrel()
    #    if kret == 0:
    #        Path('{0}/kestrel.rt'.format(completion_path)).touch()
    #if kret != 0:
    #    raise RuntimeError('Kestrel failed to complete; Exiting MARs')
    #else:
    #    main_logger.debug('Kestrel stats completed')

    #Filer  and annotate variant calls
    main_logger.debug('Annotating variants')
    annotate = Vcf.Annotate()
    gvcf_path = annotate.getAnnotation(bed_path, gvcf_path, ref_path, out_path, bam_path)
    vcf_path = annotate.getAnnotation(bed_path, vcf_path, ref_path, out_path, bam_path)
    main_logger.debug('Filetering low quality variants and merging GATK and Samtools calls')
    gvcf_file = Vcf.Reader(gvcf_path)
    svcf_file = Vcf.Reader(vcf_path)
    merged_vcf = Vcf.Merge(gvcf_file, svcf_file, out_path).merge()
    config = dict()
    summary = Summary(ref_path, bed_path, voi_path, out_dir, config)
    var_sum = summary.getVarStats(merged_vcf)
    main_logger.info('Total variants : {0}; Verified calls : {1}; Exonic : {2}; Intronic : {3}; Synonymous : {4}; Non Synonymous : {5}; Transition : {6}; Transversion : {7}'.format(
                        var_sum[0], var_sum[1], var_sum[2], var_sum[3], var_sum[4], var_sum[5], var_sum[6], var_sum[7]))
    return(merged_vcf, 0)

def marsBatch(bbduk_path, aligner_path, smt_path, bft_path, gatk_path,
              inp_path, ref_path, adp_path, bed_path, out_dir, aligner,
              pic_path, voi_path, java_path, sra_path, verbose):
    #Creating logger for nest
    logger = logging.getLogger('NeST')
    logger.setLevel(logging.DEBUG)
    #Create output paths for the run
    if not os.path.exists(os.path.abspath(out_dir)):
        os.mkdir(os.path.abspath(out_dir))
    # Creating a file handler which logs even debug messages
    fh = logging.FileHandler('{0}/nest.log'.format(os.path.abspath(out_dir)))
    if verbose:
        fh.setLevel(logging.DEBUG)
    else:
        fh.setLevel(logging.INFO)
    # Creating a console handler to log info messages
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    # Create formatter and add it to the handlers
    formatter = logging.Formatter('{asctime} - {name} - {levelname} - {message}', style="{")
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    # Add the handlers to the logger
    logger.addHandler(fh)
    logger.addHandler(ch)
    #Create file and console handlers for MaRS
    logger.info('Gathering input information from input path.')
    prep = Prepper(inp_path, sra_path)
    config = prep.prepInputs()
    logger.info('Running MaRS on {0} experiments'.format(len(config)))
    #summary = Summary(ref_path, bed_path, voi_path, out_dir)
    samples = config.keys()
    pools = Pool(5)
    rone_list = list()
    rtwo_list = list()
    name_list = list()
    for samples in config:
        name_list.append(config[samples].sample)
        rone_list.append(config[samples].files[0])
        rtwo_list.append(config[samples].files[1])


    vcf_list = pools.map(main, zip(repeat(bbduk_path), repeat(aligner_path),
                repeat(smt_path), repeat(bft_path), repeat(gatk_path),
                rone_list, rtwo_list, repeat(ref_path), repeat(adp_path),
                repeat(bed_path), repeat(out_dir), repeat(aligner),
                repeat(pic_path), name_list, repeat(voi_path),
                repeat(java_path)))
    
    if voi_path is not None:
        logger.info('Summarizing variant calls from all {0} experiments'.format(len(config)))
        summary = Summary(ref_path, bed_path, voi_path, out_dir, config)
        #Sumarize variants of intrest
        summary.getSummary()
    elif voi_path is None:
        logging.info('Variant of interest file not provided, skipping Summarize')
    return(0)

if __name__ == '__main__':
    #Define deffault paths and aligner informations
    def_path = "{0}/lib".format(os.path.abspath(os.path.dirname(os.path.realpath(__file__))))
    ref_def_path = "{0}/ref".format(os.path.abspath(os.path.dirname(os.path.realpath(__file__))))
    bbduk_def = 'bbduk.sh' #"{0}/bbmap/bbduk.sh".format(def_path)
    bbmap_def = 'bbmap.sh' #"{0}/bbmap/bbmap.sh".format(def_path)
    bwa_def = 'bwa' #"{0}/bwa/bwa".format(def_path)
    bowtie_def = 'bowtie2' #"{0}/bowtie2/bowtie2".format(def_path)
    snap_def = 'snap-alinger' #"{0}/snap/snap-aligner".format(def_path)
    smt_def = 'samtools' #"{0}/samtools/samtools".format(def_path)
    bft_def = 'bcftools' #"{0}/bcftools/bcftools".format(def_path)
    gatk_def = 'gatk' #"{0}/GenomeAnalysisTK.jar".format(def_path)
    pic_def = 'picard' #"{0}/picard.jar".format(def_path)
    sra_def = 'fastq-dump' #'{0}/sratoolkit/bin/fastq-dump'.format(def_path)
    voi_def = None #'{0}/Reportable_SNPs.csv'.format(ref_def_path)
    #if 'java version "1.8.' in str(subprocess.check_output(["java", "-version"], stderr=subprocess.STDOUT).decode('UTF-8').split('\n')[0]):
    java_def = 'java'
    #else:
    #    java_def = "{0}/jdk/bin/java".format(def_path)
    aligner_def = {'bwa' : bwa_def, 'snap' : snap_def, 'bowtie2': bowtie_def, 'bbmap': bbmap_def}
    #Get arguments
    parser = argparse.ArgumentParser(prog='NeST')
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
    parser.add_argument('--varofint', dest='voi_path', type=str, default=voi_def,
                        help='Path to variant of interest')
    parser.add_argument('--verbose', action='store_true', 
                        help='Increase verbosity of log file')                        
    args = parser.parse_args()

    #Validate parsed arguments
    if args.aligner_path is None:
        args.aligner_path = aligner_def[args.aligner]

    if not os.path.exists(args.out_path):
        os.mkdir(args.out_path)

    #Check if the run command is for batch mode analysis or single sample
    #analysis.
    #If inp_path is empty and rone_path is not, then the experiment is a
    #single sample experiment.
    status = marsBatch(args.bbduk_path, args.aligner_path, args.smt_path,
                args.bft_path, args.gatk_path, args.inp_path, args.ref_path,
                args.adp_path, args.bed_path, args.out_path, args.aligner,
                args.pic_path, args.voi_path, java_def, sra_def, args.verbose)
