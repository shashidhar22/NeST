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
from pyamd.bbduk import QualCheck
from pyamd.alignment import Bwa
from pyamd.alignment import Bowtie
from pyamd.alignment import BBMap
from pyamd.alignment import Snap
from pyamd.samtools import Samtools
from pyamd.gatk import GenAnTK
from pyamd.gatk import Picard
from pyamd.kestrel import KestrelVar
#from pyamd.annotater import Annotate
from pyamd.kestrel import kes_runner
from pyamd.summarize import Summary
from pyamd.prepinputs import Prepper
from pyamd.parsers.vcf import Vcf

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
    main_logger = logging.getLogger('Kookaburra.{0}'.format(sam_name))
    #Check if files are present
    #sam_name = config[samples].sample
    #rone_path = config[samples].files[0]
    #rtwo_path = config[samples].files[1]
    out_path = '{0}/{1}'.format(os.path.abspath(out_dir), sam_name)
    if not os.path.exists(out_path):
        os.mkdir(out_path)
    #main_logger.info('Analyzing sample : {0}'.format(sam_name))


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
        bam_path = '{0}/{1}_fixmate.bam'.format(out_path, base)
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
        bam_path = '{0}/{1}_sorted.bam'.format(out_path, base)
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

    rgadder = Picard(java_path, pic_path, out_path)
    if os.path.exists('{0}/readgroup.rt'.format(completion_path)):
        base = os.path.splitext(os.path.basename(bam_path))[0]
        bam_path = '{0}/{1}_RG.bam'.format(out_path, base)
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
        bcf_path = '{0}/variants.bcf'.format(out_path)
        pret = 0
        main_logger.debug('Skipping Pileup')
    else:
        bcf_path, pret = varengine.pileup(ref_path, bam_path)
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
        vcf_path = '{0}/{1}_variants.vcf'.format(out_path, sam_name)
        bret = 0
        main_logger.debug('Skipping bcfcall')
    else:
        vcf_path, bret = varengine.bcftools(bcf_path, bed_path, sam_name)
        main_logger.debug('Running Bcftools call')
        if bret == 0:
            Path('{0}/bcfcall.rt'.format(completion_path)).touch()

    if bret != 0:
        raise RuntimeError('Bcftools call failed to complete; Exiting MARs')
    else:
        main_logger.debug('Bcftools call completed')

    #if os.path.exists('{0}/stats.rt'.format())
    #stats_path, bret = varengine.bcfstats(vcf_path, ref_path)
    #main_logger.debug('Running Bcftools stats')
    #if bret != 0:
    #    raise RuntimeError('Bcftools stats failed to complete; Exiting MARs')
    #else:
    #    main_logger.debug('Bcftools stats completed')

    #Call GATK HaplotypeCaller to generate VCF files
    varcaller = GenAnTK(gatk_path, out_path, java_path)
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
    gvcf_path = annotate.getAnnotation(bed_path, gvcf_path, ref_path, out_path)
    vcf_path = annotate.getAnnotation(bed_path, vcf_path, ref_path, out_path)
    main_logger.debug('Filetering low quality variants and merging GATK and Samtools calls')
    gvcf_file = Vcf.Reader(gvcf_path)
    svcf_file = Vcf.Reader(vcf_path)
    merge_vcf = Vcf.Merge(gvcf_file, svcf_file)
    merged_vcf = merge_vcf.merge(out_path)
#    merged_vcf = annotate.iterVcf(bed_path, merged_vcf, sam_name, ref_path, 'merged'7)
#    gatk_vcf = annotate.iterVcf(bed_path, gvcf_path, sam_name, ref_path, 'gatk')
#    samtools_vcf = annotate.iterVcf(bed_path, vcf_path , sam_name, ref_path, 'samtools')
    summary = Summary(ref_path, bed_path, voi_path, out_dir)
    var_sum = summary.getVarStats(merged_vcf)
    main_logger.info('Total variants : {0}; Verified calls : {1}; Exonic : {2}; Intronic : {3}; Synonymous : {4}; Non Synonymous : {5}; Transition : {6}; Transversion : {7}'.format(
                        var_sum[0], var_sum[1], var_sum[2], var_sum[3], var_sum[4], var_sum[5], var_sum[6], var_sum[7]))
    return(merged_vcf, 0)

def marsBatch(bbduk_path, aligner_path, smt_path, bft_path, gatk_path,
              inp_path, ref_path, adp_path, bed_path, out_dir, aligner,
              pic_path, voi_path, java_path, sra_path):
    #Creating logger for pyamd
    logger = logging.getLogger('Kookaburra')
    logger.setLevel(logging.DEBUG)
    #Create output paths for the run
    if not os.path.exists(os.path.abspath(out_dir)):
        os.mkdir(os.path.abspath(out_dir))
    # Creating a file handler which logs even debug messages
    fh = logging.FileHandler('{0}/kookaburra.log'.format(os.path.abspath(out_dir)))
    fh.setLevel(logging.DEBUG)
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
                repeat(pic_path), name_list, repeat(voi_path),
                repeat(java_path)))

    logger.info('Summarizing variant calls from all {0} experiments'.format(len(config)))
    summary = Summary(ref_path, bed_path, voi_path, out_dir)
    #Sumarize variants of intrest
    exp_voi = summary.getRepSnps()
    exp_voi = summary.getDepthStats(exp_voi)
    exp_voi = exp_voi.reset_index(level=1)
    #exp_voi.drop_duplicates(subset='Variant', inplace=True)
    exp_voi[['Gene_name', 'RefAA_sym', 'AAPos_sort', 'AltAA_sym']] = exp_voi['Variant'].str.extract('(?P<Gene_name>[a-zA-Z0-9]+):(?P<RefAA_sym>[a-zA-Z]?)(?P<AAPos_sort>[0-9]+)(?P<AltAA_sym>[a-zA-Z]?)', expand=True)
    exp_voi['Sample_name'] = exp_voi.index
    exp_voi['AAPos_sort'] = pd.to_numeric(exp_voi['AAPos_sort'])
    exp_voi.sort_values(['Sample_name', 'Gene_name', 'AAPos_sort'], inplace=True)
    exp_voi.drop(labels=['Sample_name', 'Gene_name', 'RefAA_sym', 'AAPos_sort',
                  'AltAA_sym','Sample'], axis=1, inplace=True)
    exp_voi.to_csv('{0}/Study_variants.csv'.format(out_dir))

    exp_af = exp_voi.pivot(exp_voi.index, 'Variant')['AF'].transpose()
    exp_af['Variant'] = exp_af.index
    exp_af[['Gene_name', 'RefAA_sym', 'AAPos_sort', 'AltAA_sym']] = exp_af['Variant'].str.extract('(?P<Gene_name>[a-zA-Z0-9]+):(?P<RefAA_sym>[a-zA-Z]?)(?P<AAPos_sort>[0-9]+)(?P<AltAA_sym>[a-zA-Z]?)', expand=True)
    exp_af['AAPos_sort'] = pd.to_numeric(exp_af['AAPos_sort'])
    exp_af.sort_values(['Gene_name', 'AAPos_sort'], inplace=True)
    exp_af.drop(labels=['Variant', 'Gene_name', 'RefAA_sym', 'AAPos_sort',
                  'AltAA_sym'], axis=1, inplace=True)
    af_mask = exp_af.isnull()
    exp_af.to_csv('{0}/Study_reportable_variants_allele_frequency.csv'.format(out_dir))
#    summary.plotHeatMap(exp_af, 'voi_af', af_mask)
    exp_dp = exp_voi.pivot(exp_voi.index, 'Variant')['DP'].transpose()
    exp_dp['Variant'] = exp_dp.index
    exp_dp[['Gene_name', 'RefAA_sym', 'AAPos_sort', 'AltAA_sym']] = exp_dp['Variant'].str.extract('(?P<Gene_name>[a-zA-Z0-9]+):(?P<RefAA_sym>[a-zA-Z]?)(?P<AAPos_sort>[0-9]+)(?P<AltAA_sym>[a-zA-Z]?)', expand=True)
    exp_dp['AAPos_sort'] = pd.to_numeric(exp_dp['AAPos_sort'])
    exp_dp.sort_values(['Gene_name', 'AAPos_sort'], inplace=True)
    exp_dp.drop(labels=['Variant', 'Gene_name', 'RefAA_sym', 'AAPos_sort',
                  'AltAA_sym'], axis=1, inplace=True)
    dp_mask = exp_dp.isnull()
    exp_dp.to_csv('{0}/Study_reportable_variants_depth.csv'.format(out_dir))
#    summary.plotHeatMap(exp_dp, 'voi_dp', dp_mask)
#    summary.plotCountPlot(exp_af, 'voi')
    #Summarize novel variants
    exp_nov = summary.getNovSnps()
    exp_nov = summary.getNovDepthStats(exp_nov)
    exp_nov = exp_nov.reset_index(level=1)
    exp_nov[['Gene_name', 'RefAA_sym', 'AAPos_sort', 'AltAA_sym']] = exp_nov['Variant'].str.extract('(?P<Gene_name>[a-zA-Z0-9]+):(?P<RefAA_sym>[a-zA-Z]?)(?P<AAPos_sort>[0-9]+)(?P<AltAA_sym>[a-zA-Z]?)', expand=True)
    exp_nov['Sample_name'] = exp_nov.index
    exp_nov['AAPos_sort'] = pd.to_numeric(exp_nov['AAPos_sort'])
    exp_nov.sort_values(['Sample_name', 'Gene_name', 'AAPos_sort'], inplace=True)
    exp_nov.drop(labels=['Sample_name', 'Gene_name', 'RefAA_sym', 'AAPos_sort',
                  'AltAA_sym', 'Sample'], axis=1, inplace=True)
    exp_nov.to_csv('{0}/Study_novel_exonic_variants.csv'.format(out_dir))
    #Separate and capture Intron and exonic variants
    exp_nov_af = exp_nov.loc[:,['Variant', 'AF']]
    exp_nov_af[['Gene_name', 'RefAA_sym', 'AAPos_sort', 'AltAA_sym']] = exp_nov_af['Variant'].str.extract('(?P<Gene_name>[a-zA-Z0-9]+):(?P<RefAA_sym>[a-zA-Z]?)(?P<AAPos_sort>[0-9]+)(?P<AltAA_sym>[a-zA-Z]?)', expand=True)
    exp_nov_af['AAPos_sort'] = pd.to_numeric(exp_nov_af['AAPos_sort'])
    exp_nov_af.sort_values(['Gene_name', 'AAPos_sort'], inplace=True)
    exp_nov_af.drop(labels=['Variant', 'Gene_name', 'RefAA_sym', 'AAPos_sort',
                  'AltAA_sym'], axis=1, inplace=True)
    exp_nov_af.to_csv('{0}/Study_novel_variants_alele_frequency.csv'.format(out_dir))
    exp_nov_dp = exp_nov.loc[:,['Variant', 'DP']]
    exp_nov_dp['Variant'] = exp_nov_dp.index
    exp_nov_dp[['Gene_name', 'RefAA_sym', 'AAPos_sort', 'AltAA_sym']] = exp_nov_dp['Variant'].str.extract('(?P<Gene_name>[a-zA-Z0-9]+):(?P<RefAA_sym>[a-zA-Z]?)(?P<AAPos_sort>[0-9]+)(?P<AltAA_sym>[a-zA-Z]?)', expand=True)
    exp_nov_dp['AAPos_sort'] = pd.to_numeric(exp_nov_dp['AAPos_sort'])
    exp_nov_dp.sort_values(['Gene_name', 'AAPos_sort'], inplace=True)
    exp_nov_dp.drop(labels=['Variant', 'Gene_name', 'RefAA_sym', 'AAPos_sort',
                  'AltAA_sym'], axis=1, inplace=True)
    exp_nov_dp.to_csv('{0}/Study_novel_variants_depth.csv'.format(out_dir))
    exp_intron = summary.getIntronTables()
    exp_intron = exp_intron.reset_index()

    #print(exp_intron.index)
    #exp_intron.reset_index(level=1)
#    print(exp_intron.head())
    exp_intron[['Gene_name', 'RefAA_sym', 'AAPos_sort', 'AltAA_sym']] = exp_intron['Variant'].str.extract('(?P<Gene_name>[a-zA-Z0-9]+):(?P<RefAA_sym>[a-zA-Z]?)(?P<AAPos_sort>[0-9]+)(?P<AltAA_sym>[a-zA-Z]?)', expand=True)
    #exp_intron['']
    exp_intron['AAPos_sort'] = pd.to_numeric(exp_intron['AAPos_sort'])
    exp_intron.sort_values(['Sample', 'Gene_name', 'AAPos_sort'], inplace=True)
    exp_intron.drop(labels=['Gene_name', 'RefAA_sym', 'AAPos_sort',
                  'AltAA_sym' ], axis=1, inplace=True)
    exp_intron.sort_index().reset_index(drop=True).to_csv('{0}/Study_novel_intronic_variants.csv'.format(out_dir), index=False)
    # Plot using Rscript
    logger.info('Plotting Depth Per SNP')
    dcmd = ['Rscript', 'pyamd/Rscripts/DepthPerReportSNP.R', '-i',
            '{0}/Study_reportable_variants_depth.csv'.format(out_dir), '-o',
            '{0}/Study_depth.pdf'.format(out_dir)]
    drun = subprocess.Popen(dcmd, shell=False,
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    drun.wait()
    if drun.returncode != 0:
        logger.error('Failed to execute DepthPerReportSNP.R')
        logger.error(' '.join(dcmd))

    logger.info('Plotting Reportable SNPs Frequency')
    acmd = ['Rscript', 'pyamd/Rscripts/reportableSNPsFreq.R', '-i',
            'Study_variants.csv'.format(out_dir), '-r',
            'ref/Reportable_SNPs.csv', '-o', '{0}/'.format(out_dir)]
    arun = subprocess.Popen(acmd, shell=False,
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    arun.wait()
    if arun.returncode != 0:
        logger.error('Failed to execute reportableSNPsFreq.R')
        logger.error(' '.join(acmd))
    logger.info('Plotting Novel Exonic Non-Synonymous SNPs')
    nenscmd = ['Rscript', 'pyamd/Rscripts/NovelExonicNonSynSNPs.R', '-i',
            'Study_novel_exonic_variants.csv'.format(out_dir),
            '-o', '{0}/'.format(out_dir)]
    nensrun = subprocess.Popen(nenscmd, shell=False,
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    nensrun.wait()
    if nensrun.returncode != 0:
        logger.error('Failed to execute NovelExonicNonSynSNPs.R')
        logger.error(' '.join(nenscmd))

    logger.info('Plotting Novel Exonic Synonymous SNPs')
    nescmd = ['Rscript', 'pyamd/Rscripts/NovelExonicSynSNPs.R', '-i',
            'Study_novel_exonic_variants.csv'.format(out_dir),
            '-o', '{0}/'.format(out_dir)]
    nesrun = subprocess.Popen(nescmd, shell=False,
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    nesrun.wait()
    if nesrun.returncode != 0:
        logger.error('Failed to execute NovelExonicSynSNPs.R')
        logger.error(' '.join(acmd))

    logger.info('Plotting Novel Intronic SNPs')
    nicmd = ['Rscript', 'pyamd/Rscripts/NovelIntronicSNPs.R', '-i',
            'Study_novel_intronic_variants.csv'.format(out_dir),
            '-o', '{0}/'.format(out_dir)]
    nirun = subprocess.Popen(nicmd, shell=False,
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    nirun.wait()
    if nirun.returncode != 0:
        logger.error('Failed to execute NovelIntronicSNPs.R')
        logger.error(' '.join(nicmd))

    #os.remove('{0}/Reportable_SNPs_Report.csv'.format(out_dir))
    os.remove('{0}/novel_SNPs_exonic_syn.csv'.format(out_dir))
    os.remove('{0}/novel_SNPs_intronic.csv'.format(out_dir))
    os.remove('{0}/novel_SNPs_exonic_nonsyn.csv'.format(out_dir))
    os.remove('{0}/Study_novel_exonic_variants_filtered.csv'.format(out_dir))
    os.remove('{0}/Study_novel_intronic_variants_filtered.csv'.format(out_dir))
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
    voi_def = '{0}/Reportable_SNPs.csv'.format(ref_def_path)
    #if 'java version "1.8.' in str(subprocess.check_output(["java", "-version"], stderr=subprocess.STDOUT).decode('UTF-8').split('\n')[0]):
    java_def = 'java'
    #else:
    #    java_def = "{0}/jdk/bin/java".format(def_path)
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
    parser.add_argument('--varofint', dest='voi_path', type=str, default=voi_def,
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
    status = marsBatch(args.bbduk_path, args.aligner_path, args.smt_path,
                args.bft_path, args.gatk_path, args.inp_path, args.ref_path,
                args.adp_path, args.bed_path, args.out_path, args.aligner,
                args.pic_path, args.voi_path, java_def, sra_def)
