#! /usr/bin/python3
import os
import sys
import glob
import time
import shutil
import logging
import argparse
import subprocess
from itertools import repeat
from multiprocessing import Pool
from pyamd.bbduk import QualCheck
from pyamd.alignment import Bwa
from pyamd.alignment import Bowtie
from pyamd.alignment import BBMap
from pyamd.alignment import Snap
from pyamd.samtools import Samtools
from pyamd.reader import Reader
from pyamd.gatk import GenAnTK
from pyamd.gatk import Picard
from pyamd.annotater import Annotate
from pyamd.kestrel import kes_runner
from pyamd.filter import filterer
from pyamd.summarize import Summary
from pyamd.prepinputs import Prepper

loggermain = logging.getLogger(__name__)
loggermain.setLevel(logging.INFO)
chmain = logging.StreamHandler()
chmain.setLevel(logging.ERROR)
formattermain = logging.Formatter('%(asctime)s:%(name)s:%(levelname)s:%(message)s')
chmain.setFormatter(formattermain)
loggermain.addHandler(chmain)

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
    kes_path = arguments[12]
    kan_path = arguments[13]
    pic_path = arguments[14]
    sam_name = arguments[15]
    voi_path = arguments[16]
    java_path = arguments[17]
    #Setup logging
    #Check if files are present
    #sam_name = config[samples].sample
    #rone_path = config[samples].files[0]
    #rtwo_path = config[samples].files[1]
    out_path = '{0}/{1}'.format(os.path.abspath(out_dir), sam_name)
    if not os.path.exists(out_path):
        os.mkdir(out_path)
    #loggermain.info('Analyzing sample : {0}'.format(sam_name))


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
    loggermain.debug('Running BBDuk')
    bbduk = QualCheck(bbduk_path, adp_path, out_path, java_path)
    rone_path, rtwo_path, bret = bbduk.bbduk(rone_path, rtwo_path)
    if bret != 0:
        raise RuntimeError('BBDuk failed to complete; Exiting MARs')
    else:
        loggermain.debug('BBDuk completed')

    if aligner == 'bwa':
        #Call BWA
        loggermain.debug('Running BWA')
        bwa = Bwa(alinger_path, out_path, ref_path)
        sam_path, mret = bwa.bwamem(rone_path, rtwo_path)
        if mret != 0:
            raise RuntimeError('Bwa mem failed to complete; Exiting MARs')
        else:
            loggermain.debug('BWA completed')

    elif aligner == 'bowtie2':
        #Call Bowtie2
        loggermain.debug('Running Bowtie2')
        bowtie = Bowtie(alinger_path, out_path, ref_path)
        sam_path, mret = bowtie.bowtie(rone_path, rtwo_path)
        if mret != 0:
            raise RuntimeError('Bowtie2 failed to complete; Exiting MARs')
        else:
            loggermain.debug('Bowtie2 completed')

    elif aligner == 'snap':
        #Call Snap
        loggermain.debug('Running Snap')
        snap = Snap(alinger_path, out_path, ref_path)
        sam_path, mret = snap.snap(rone_path, rtwo_path)
        if mret != 0:
            raise RuntimeError('Snap failed to complete; Exiting MARs')
        else:
            loggermain.debug('Snap completed')

    elif aligner == 'bbmap':
        #Call Bbmap
        loggermain.debug('Running BBMap')
        bbmap = BBMap(alinger_path, out_path, ref_path)
        sam_path, mret = bbmap.bbmap(rone_path, rtwo_path)
        if mret != 0:
            raise RuntimeError('BBMap failed to complete; Exitinign MARs')
        else:
            loggermain.debug('BBMap completed')


    #Fix mate information, sort files and add read groups
    varengine = Samtools(smt_path, bft_path, out_path)
    bam_path, fret = varengine.fixmate(sam_path)
    loggermain.debug('Running Samtools fixmate')
    if fret != 0:
        raise RuntimeError('Samtools fixmate failed to complete; Exiting MARs')
    else:
        loggermain.debug('Samtools fixmate completed')

    bam_path, sret = varengine.sort(sam_path)
    loggermain.debug('Running Samtools sort')
    if sret != 0:
        raise RuntimeError('Samtools sort failed to complete; Exiting MARs')
    else:
        loggermain.debug('Samtools sort completed')

    rgadder = Picard(java_path, pic_path, out_path)
    bam_path, aret = rgadder.picard(bam_path, sam_name)
    loggermain.debug('Running Picard AddOrReplaceReadGroups')
    if aret != 0:
        raise RuntimeError('Picard AddOrReplaceReadGroups failed to complete; Exiting MARs')
    else:
        loggermain.debug('Picard AddOrReplaceReadGroups completed')

    #Run samtools mpileup, bcftools index, call and stats to generate VCF files
    bcf_path, pret = varengine.pileup(ref_path, bam_path)
    loggermain.debug('Running Samtools mpileup')
    if pret != 0:
        raise RuntimeError('Samtools mpileup failed to complete; Exiting MARs')
    else:
        loggermain.debug('Samtools mpileup completed')

    bret = varengine.bcfindex(bcf_path)
    loggermain.debug('Running Bcftools index')
    if bret != 0:
        raise RuntimeError('Bcftools index failed to complete; Exiting MARs')
    else:
        loggermain.debug('Bcftools index completed')

    vcf_path, bret = varengine.bcftools(bcf_path, bed_path, sam_name)
    loggermain.debug('Running Bcftools call')
    if bret != 0:
        raise RuntimeError('Bcftools call failed to complete; Exiting MARs')
    else:
        loggermain.debug('Bcftools call completed')

    stats_path, bret = varengine.bcfstats(vcf_path, ref_path)
    loggermain.debug('Running Bcftools stats')
    if bret != 0:
        raise RuntimeError('Bcftools stats failed to complete; Exiting MARs')
    else:
        loggermain.debug('Bcftools stats completed')

    #Call GATK HaplotypeCaller to generate VCF files
    varcaller = GenAnTK(gatk_path, out_path, java_path)
    loggermain.debug('Running GATK HaplotypeCaller')
    gvcf_path, gret = varcaller.hapCaller(bam_path, ref_path, sam_name)
    if gret != 0:
        raise RuntimeError('GATK HaplotypeCaller failed to complete; Exiting MARs')
    else:
        loggermain.debug('GATK HaplotypeCaller stats completed')

    #Filer  and annotate variant calls
    loggermain.debug('Filetering low quality variants and merging GATK and Samtools calls')
    merged_vcf = filterer(gvcf_path, vcf_path, sam_name, out_path)
    loggermain.debug('Annotating variants')
    annotate = Annotate(out_path)
    merged_vcf = annotate.iterVcf(bed_path, merged_vcf, sam_name, ref_path, 'merged')
    gatk_vcf = annotate.iterVcf(bed_path, gvcf_path, sam_name, ref_path, 'gatk')
    samtools_vcf = annotate.iterVcf(bed_path, vcf_path , sam_name, ref_path, 'samtools')
    summary = Summary(ref_path, bed_path, voi_path, out_dir)
    var_sum = summary.getVarStats(merged_vcf)
    loggermain.info('Finished analyzing sample : {8} \n Total variants : {0}; Verified calls : {1}; Exonic : {2}; Intronic : {3}; Synonymous : {4}; Non Synonymous : {5}; Transition : {6}; Transversion : {7}'.format(
                        var_sum[0], var_sum[1], var_sum[2], var_sum[3], var_sum[4], var_sum[5], var_sum[6], var_sum[7], sam_name))
    return(merged_vcf, 0)

def marsBatch(bbduk_path, aligner_path, smt_path, bft_path, gatk_path,
              inp_path, ref_path, adp_path, bed_path,
              out_dir, aligner, kes_path, kan_path, pic_path, voi_path, java_path):
    #Create loggermain for MaRS
    #Create output paths for the run
    if not os.path.exists(os.path.abspath(out_dir)):
        os.mkdir(os.path.abspath(out_dir))
    #Create file and console handlers for MaRS
    loggermain.info('Gathering input information from input path.')
    prep = Prepper(inp_path)
    config = prep.prepInputs()
    loggermain.info('Running MaRS on {0} experiments'.format(len(config)))
    #summary = Summary(ref_path, bed_path, voi_path, out_dir)
    samples = config.keys()
    pools = Pool(4)
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
                repeat(kes_path), repeat(kan_path), repeat(pic_path), name_list,
                repeat(voi_path), repeat(java_path)))

    loggermain.info('Summarizing variant calls from all {0} experiments'.format(len(config)))
    summary = Summary(ref_path, bed_path, voi_path, out_dir)
    #Sumarize variants of intrest
    exp_voi = summary.getRepSnps()
    exp_voi = summary.getDepthStats(exp_voi)
    exp_voi = exp_voi.reset_index(level=1)
    #exp_voi.drop_duplicates(subset='Variant', inplace=True)
    exp_voi.sort_index().to_excel('{0}/Study_variants.xlsx'.format(out_dir))

    exp_af = exp_voi.pivot(exp_voi.index, 'Variant')['AF'].transpose()
    af_mask = exp_af.isnull()
    exp_af.sort_index().to_excel('{0}/Study_al_freq.xlsx'.format(out_dir))
    summary.plotHeatMap(exp_af, 'voi_af', af_mask)
    exp_dp = exp_voi.pivot(exp_voi.index, 'Variant')['DP'].transpose()
    dp_mask = exp_dp.isnull()
    exp_dp.sort_index().to_excel('{0}/Study_depth.xlsx'.format(out_dir))
    summary.plotHeatMap(exp_dp, 'voi_dp', dp_mask)
    summary.plotCountPlot(exp_af, 'voi')
    #Summarize novel variants
    exp_nov = summary.getNovSnps()
    exp_nov = summary.getNovDepthStats(exp_nov)
    exp_nov = exp_nov.reset_index(level=1)
    #Separate and capture Intron and exonic variants
    exp_nov_af = exp_nov.loc[:,['Variant', 'AF']]
    exp_nov_af.sort_index().to_excel('{0}/Study_novel_var_af.xlsx'.format(out_dir))
    exp_nov_dp = exp_nov.loc[:,['Variant', 'DP']]
    exp_nov_dp.sort_index().to_excel('{0}/Study_novel_var_dp.xlsx'.format(out_dir))
    exp_nov.sort_index().to_excel('{0}/Study_novel_exonic_variants.xlsx'.format(out_dir))
    exp_intron = summary.getIntronTables()
    exp_intron.sort_index().to_excel('{0}/Study_novel_intronic_variants.xlsx'.format(out_dir))
    return(0)

if __name__ == '__main__':
    #Define deffault paths and aligner informations
    def_path = "{0}/lib".format(os.path.abspath(os.path.dirname(os.path.realpath(__file__))))
    bbduk_def = "{0}/bbmap/bbduk.sh".format(def_path)
    bbmap_def = "{0}/bbmap/bbmap.sh".format(def_path)
    bwa_def = "{0}/bwa-0.7.12/bwa".format(def_path)
    bowtie_def = "{0}/bowtie2-2.3.0/bowtie2".format(def_path)
    snap_def = "{0}/snap/snap-aligner".format(def_path)
    smt_def = "{0}/samtools-1.3.1/samtools".format(def_path)
    bft_def = "{0}/bcftools-1.3.1/bcftools".format(def_path)
    gatk_def = "{0}/GenomeAnalysisTK.jar".format(def_path)
    pic_def = "{0}/picard.jar".format(def_path)
    if 'java version "1.8.' in str(subprocess.check_output(["java", "-version"], stderr=subprocess.STDOUT).decode('UTF-8').split('\n')[0]):
        java_def = 'java'
    else:
        java_def = "{0}/jdk1.8.0_131/bin/java".format(def_path)
    aligner_def = {'bwa' : bwa_def, 'snap' : snap_def, 'bowtie2': bowtie_def, 'bbmap': bbmap_def}
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

    #Validate parsed arguments
    if args.aligner_path is None:
        args.aligner_path = aligner_def[args.aligner]

    if not os.path.exists(args.out_path):
        os.mkdir(args.out_path)

    #Check if the run command is for batch mode analysis or single sample
    #analysis.
    #If inp_path is empty and rone_path is not, then the experiment is a
    #single sample experiment.
    if args.inp_path == None and args.rone_path != None:
        if args.sam_name == None:
            sam_name = os.path.splitext(os.path.basename(args.rone_path))[0]
        else:
            sam_name = args.sam_name
        status = main(args.bbduk_path, args.aligner_path, args.smt_path,
                    args.bft_path, args.gatk_path, args.rone_path,
                    args.rtwo_path, args.ref_path, args.adp_path, args.bed_path,
                    args.out_path, args.aligner, args.kes_path, args.kan_path,
                    args.pic_path, sam_name, java_def)
    elif args.inp_path != None and args.rone_path == None:
        status = marsBatch(args.bbduk_path, args.aligner_path, args.smt_path,
                    args.bft_path, args.gatk_path, args.inp_path, args.ref_path,
                    args.adp_path, args.bed_path, args.out_path, args.aligner,
                    args.kes_path, args.kan_path, args.pic_path, args.voi_path,
                    java_def)
