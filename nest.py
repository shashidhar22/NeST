#! /usr/bin/python3
#####comment: written in OOP and classes
#####comment: import tools
####reference: https://www.webopedia.com/TERM/O/object_oriented_programming_OOP.html
####comment: OOP: Object-oriented programming (OOP) refers to a type of computer 
####programming (software design) in which programmers define the data type of a data structure, 
####and also the types of operations (functions) that can be applied to the data structure.
####In this way, the data structure becomes an object that includes both data and functions. 
####In addition, programmers can create relationships between one object and another. For example, 
####objects can inherit characteristics from other objects.
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
#####comment: import bioinformatics tools
from nest.alignment import Bwa #####comment: import Bwa aligner
from nest.alignment import Bowtie #####comment: import Bowtie aligner
from nest.alignment import BBMap #####comment: import BBMAP
from nest.alignment import Snap
from nest.samtools import Samtools
from nest.gatk import GenAnTK
from nest.gatk import Picard
from nest.gatk import FreeBayes
from nest.kestrel import KestrelVar
#from nest.annotater import Annotate
from nest.kestrel import kes_runner
from nest.summarize import Summary
from nest.prepinputs import Prepper
from nest.parsers.vcfReader import Reader 
from nest.parsers.vcfmerge import Merge
from nest.parsers.vcfannotate import Annotate
from nest.parsers.vcfwriter import Writer

#####comment: setting path for differnt moduels
#####comment: path for different tools
def main(arguments):
    bbduk_path = arguments[0]
    alinger_path = arguments[1]
    smt_path = arguments[2]
    bft_path = arguments[3]
    gatk_path = arguments[4]
    sam_name = arguments[5]
    file_list = arguments[6]
    ref_path = arguments[7]
    adp_path = arguments[8]
    bed_path = arguments[9]
    out_dir = arguments[10]
    aligner = arguments[11]
    pic_path = arguments[12]
    voi_path = arguments[13]
    java_path = arguments[14]
    sra_path = arguments[15]
    purge = arguments[16]    
    sra_list = arguments[17]
    #Setup logging
    #####comment:By logging useful data from the right places, 
    ####you can not only debug errors easily but also use the data to analyze 
    ####the performance of the application to plan for scaling or look at usage 
    ####patterns to plan for marketing.
    #Get logger for main method
    main_logger = logging.getLogger('NeST.{0}'.format(sam_name))
    main_logger.debug('Starting analysis for {0}'.format(sam_name))
    #####comment:Check if files are present
    out_path = '{0}/{1}'.format(os.path.abspath(out_dir), sam_name)
    if not os.path.exists(out_path):
        os.mkdir(out_path)

    #####comment:check paths for fastq file
    fastq_path = '{0}/RawFastq'.format(out_path)
    if not os.path.exists(fastq_path):
        os.mkdir(fastq_path)
    #####comment:Get FASTQs
    prepper =  Prepper(fastq_path, out_dir, sra_path)
    fastq_path = prepper.sra(sam_name, sra_list, file_list)
    ##Note: Generalize this, right now it will only work with SRA. This is a fix for NEJM
    rone_path = file_list[0]
    rtwo_path = file_list[1]

    #####comment:List errors when there is no specific directory input
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
    #####referemce:https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/
    #####comment:“Duk” stands for Decontamination Using Kmers. 
    #####BBDuk was developed to combine most common data-quality-related 
    #####trimming, filtering, and masking operations into a single high-performance tool. 
    #####It is capable of quality-trimming and filtering, adapter-trimming, contaminant-filtering 
    #####via kmer matching, sequence masking, GC-filtering, length filtering, entropy-filtering, 
    #####format conversion, histogram generation, subsampling, quality-score recalibration, kmer 
    #####cardinality estimation, and various other operations in a single pass. 
    main_logger.debug('Running BBDuk')
    if os.path.exists('{0}/bbduk.rt'.format(completion_path)):
        ###splittext method breaks the Text node into two nodes at the specified offset, keeping both nodes in the tree as siblings.
        ###After the split, the current node contains all the content up to the specified offset point, and a newly created node of the same type contains the remaining text. The newly created node is returned to the caller. If the original node had a parent, the new node is inserted as the next sibling of the original node. If the offset is equal to the length of the original node, the newly created node has no data.
        ###reference: https://docs.python.org/2/library/os.path.html
        ###comment: Split the pathname path into a pair (root, ext) such that root + ext == path, 
        ###and ext is empty or begins with a period and contains at most one period. 
        ###Leading periods on the basename are ignored; splitext('.cshrc') returns ('.cshrc', '').
        ###comment: before cleaning R1 and R2 of the sample
        brone = os.path.splitext(os.path.basename(rone_path))[0]
        brtwo = os.path.splitext(os.path.basename(rtwo_path))[0]
        ###comment: after cleaning R1 and R2 of the sample
        rone_path = '{0}/{1}/{2}_cleaned.fq'.format(out_path, 'CleanedFastq', brone)
        rtwo_path = '{0}/{1}/{2}_cleaned.fq'.format(out_path, 'CleanedFastq', brtwo)
        main_logger.debug('Skipping BBDuk')
        bret = 0
    else:
        ###comment: '''QualCheck class is written to filter reads from a sample, based on
        ###adapter contamination and low quality reads using BBDuk. The reads are
        ###trimmed from both sides, reads are scanned across using kmers of max length
        ###27, and minimum length 4. hdist specfies a hamming distance of 1, that is,
        ###a max of one mismatch between the kmer and the adapter sequences.
        ###Sections with average quality less 30 are trimmed. Reads smaller than 50
        ###are excluded from the study.
        ###Basically the purpose of qualcheck is to filter reads and trim it. 
        bbduk = QualCheck(bbduk_path, adp_path, out_path, java_path)
        rone_path, rtwo_path, bret = bbduk.bbduk(rone_path, rtwo_path)
        if bret == 0:
            Path('{0}/bbduk.rt'.format(completion_path)).touch()
    ###comment: When there is no path exit the program
    if bret != 0:
        raise RuntimeError('BBDuk failed to complete; Exiting MARs')
    else:
        main_logger.debug('BBDuk completed')

    ####reference: http://bio-bwa.sourceforge.net/
    ####BWA is a software package for mapping low-divergent sequences against a large reference genome, 
    ####such as the human genome. It consists of three algorithms: BWA-backtrack, BWA-SW and BWA-MEM. 
    ####The first algorithm is designed for Illumina sequence reads up to 100bp, while the rest two for 
    ####longer sequences ranged from 70bp to 1Mbp. BWA-MEM and BWA-SW share similar features such as 
    ####long-read support and split alignment, but BWA-MEM, which is the latest, is generally recommended 
    ####for high-quality queries as it is faster and more accurate. BWA-MEM also has better performance 
    ####than BWA-backtrack for 70-100bp Illumina reads. 
    if aligner == 'bwa':
        #Call BWA class 
        ###''Bwa class runs BWA mem in standard settings to align the samples reads
        ###against the reference genome
        main_logger.debug('Running BWA')
        if os.path.exists('{0}/align.rt'.format(completion_path)):
            sam_path = '{0}/alignments/output.sam'.format(out_path)
            mret = 0
            main_logger.debug('Skipping BWA')
        else:
            bwa = Bwa(alinger_path, out_path, ref_path)
            ###run bwamem function from bwa class
            sam_path, mret = bwa.bwamem(rone_path, rtwo_path)
            if mret == 0:
                Path('{0}/align.rt'.format(completion_path)).touch()
        if mret != 0:
            raise RuntimeError('Bwa mem failed to complete; Exiting MARs')
        else:
            main_logger.debug('BWA completed')

    elif aligner == 'bowtie2':
        ###reference: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
        ###comment: Bowtie 2 is an ultrafast and memory-efficient tool for aligning sequencing 
        #reads to long reference sequences. It is particularly good at aligning reads of about 50 
        #up to 100s or 1,000s of characters, and particularly good at aligning to relatively long 
        #(e.g. mammalian) genomes. Bowtie 2 indexes the genome with an FM Index to keep its memory 
        #footprint small: for the human genome, its memory footprint is typically around 3.2 GB. 
        #Bowtie 2 supports gapped, local, and paired-end alignment modes. 
        #Call Bowtie2
        main_logger.debug('Running Bowtie2')
        if os.path.exists('{0}/aling.rt'.format(completion_path)):
            sam_path = '{0}/alignments/output.sam'.format(out_path)
            mret = 0
            main_logger.debug('Skipping Bowtie2')
        else:
            bowtie = Bowtie(alinger_path, out_path, ref_path)
            ###run bowtie function from bowtie class
            sam_path, mret = bowtie.bowtie(rone_path, rtwo_path)
            if mret == 0:
                Path('{0}/align.rt'.format(completion_path)).touch()
        if mret != 0:
            raise RuntimeError('Bowtie2 failed to complete; Exiting MARs')
        else:
            main_logger.debug('Bowtie2 completed')

    elif aligner == 'snap':
        ###reference: http://snap.cs.berkeley.edu/
        ###comment: SNAP is a new sequence aligner that is 3-20x faster and just as 
        #accurate as existing tools like BWA-mem, Bowtie2 and Novoalign. It runs on commodity 
        #x86 processors, and supports a rich error model that lets it cheaply match reads with 
        #more differences from the reference than other tools. This gives SNAP up to 2x lower 
        #error rates than existing tools (in some cases) and lets it match larger mutations 
        #that they may miss. SNAP also natively reads BAM, FASTQ, or gzipped FASTQ, and natively 
        #writes SAM or BAM, with built-in sorting, duplicate marking, and BAM indexing.
        #Call Snap
        main_logger.debug('Running Snap')
        ###run snap function from snap class
        snap = Snap(alinger_path, out_path, ref_path)
        sam_path, mret = snap.snap(rone_path, rtwo_path)
        if mret != 0:
            raise RuntimeError('Snap failed to complete; Exiting MARs')
        else:
            main_logger.debug('Snap completed')

    elif aligner == 'bbmap':
        #Call Bbmap
        ###refernece: https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/
        ###comment: BBMap is a splice-aware global aligner for DNA and RNA sequencing reads. 
        ###It can align reads from all major platforms – Illumina, 454, Sanger, 
        ###Ion Torrent, Pac Bio, and Nanopore. BBMap is fast and extremely accurate, 
        ###particularly with highly mutated genomes or reads with long indels, even 
        ###whole-gene deletions over 100kbp long. It has no upper limit to genome size or 
        ###number of contigs, and has been successfully used for mapping to an 85 gigabase soil 
        ###metagenome with over 200 million contigs. Additionally, the indexing phase is very fast 
        ###compared to other aligners.
        ###BBMap has a large array of options, described in its shell script. It can output many 
        ####different statistics files, such as an empirical read quality histogram, insert-size 
        ###distribution, and genome coverage, with or without generating a sam file. As a result, 
        ###it is useful in quality control of libraries and sequencing runs, or evaluating new 
        ###sequencing platforms. The derivative program BBSplit is also useful in binning or 
        ###refining metagenomic reads.
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

    ####reference: https://broadinstitute.github.io/picard/
    ####comment: Picard is a set of command line tools for manipulating high-throughput sequencing (HTS) 
    ###data and formats such as SAM/BAM/CRAM and VCF. These file formats are defined in the 
    ###Hts-specs repository. See especially the SAM specification and the VCF specification.
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
    ###refernece: https://gatk.broadinstitute.org/hc/en-us
    ###The GATK is the industry standard for identifying SNPs and indels in germline DNA and 
    #RNAseq data. Its scope is now expanding to include somatic short variant calling, and to 
    #tackle copy number (CNV) and structural variation (SV). In addition to the variant callers 
    #themselves, the GATK also includes many utilities to perform related tasks such as processing 
    #and quality control of high-throughput sequencing data, and bundles the popular Picard toolkit.
    ###These tools were primarily designed to process exomes and whole genomes generated with Illumina 
    #sequencing technology, but they can be adapted to handle a variety of other technologies and 
    #experimental designs. And although it was originally developed for human genetics, the GATK has 
    #since evolved to handle genome data from any organism, with any level of ploidy.
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

    #Call Freebayes to generate VCF files
    ###reference: https://github.com/ekg/freebayes
    ###freebayes is a Bayesian genetic variant detector designed to find small polymorphisms, 
    #specifically SNPs (single-nucleotide polymorphisms), indels (insertions and deletions), 
    #MNPs (multi-nucleotide polymorphisms), and complex events (composite insertion and 
    #substitution events) smaller than the length of a short-read sequencing alignment.
    ###freebayes is haplotype-based, in the sense that it calls variants based on the literal 
    #sequences of reads aligned to a particular target, not their precise alignment. This model 
    #is a straightforward generalization of previous ones (e.g. PolyBayes, samtools, GATK) which 
    #detect or report variants based on alignments. This method avoids one of the core problems 
    #with alignment-based variant detection--- that identical sequences may have multiple possible 
    #alignments:
    varcaller = FreeBayes('freebayes', out_path)
    main_logger.debug('Running Freebayes')
    if os.path.exists('{0}/freebayes.rt'.format(completion_path)):
        fvcf_path = '{0}/{1}_variants_freebayes.vcf'.format(out_path, sam_name)
        fret = 0
        main_logger.debug('Skipping Freebayes')
    else:
        fvcf_path, fret = varcaller.freeBayes(bam_path, ref_path, sam_name)
        if fret == 0:
            Path('{0}/freebayes.rt'.format(completion_path)).touch()
    if fret != 0:
        raise RuntimeError('Freebayes failed to complete; Exiting MARs')
    else:
        main_logger.debug('Freebayes stats completed')


    #Filer  and annotate variant calls
    main_logger.debug('Annotating variants')
    annotate = Annotate()
    gvcf_path = annotate.getAnnotation(bed_path, gvcf_path, ref_path, out_path, bam_path)
    vcf_path = annotate.getAnnotation(bed_path, vcf_path, ref_path, out_path, bam_path)
    fvcf_path = annotate.getAnnotation(bed_path, fvcf_path, ref_path, out_path, bam_path)
    vcf_dict = {gvcf_path: 'GATK', vcf_path: 'Samtools', fvcf_path: 'Freebayes'}
    merger = Merge(out_path, vcf_dict, ref_path)
    merged_vcf = merger.splitter(list(vcf_dict.keys()))[0]
    final_vcf= '{0}/{1}_variants_merged_annotated.vcf'.format(out_path, sam_name)
    os.rename(merged_vcf, final_vcf)
    #final_path = annotate.getAnnotation(bed_path, final_vcf, ref_path, out_path, bam_path)
    main_logger.debug('Filetering low quality variants and merging GATK and Samtools calls')
    #merged_vcf = Vcf.Merge(gvcf_file, svcf_file, out_path).merge()
    summary = Summary(ref_path, bed_path, voi_path, out_dir)
    var_sum = summary.getVarStats(final_vcf)
    main_logger.info('Total variants : {0}; Verified calls : {1}; Exonic : {2}; Intronic : {3}; Synonymous : {4}; Non Synonymous : {5}; Transition : {6}; Transversion : {7}'.format(
                        var_sum[0], var_sum[1], var_sum[2], var_sum[3], var_sum[4], var_sum[5], var_sum[6], var_sum[7]))
    if purge:
       shutil.rmtree('{0}/RawFastq'.format(out_path))
       shutil.rmtree('{0}/CleanedFastq'.format(out_path))
       alignments = glob.glob('{0}/alignments/*'.format(out_path))
       for files in alignments:
           if 'output_FM_SR_DD_RG.ba' in files:
               continue
           else:
               os.remove(files)
       vcffiles = glob.glob('{0}/*.bcf*'.format(out_path))
       for files in vcffiles:
           os.remove(files)
    return(final_vcf, 0)

def marsBatch(bbduk_path, aligner_path, smt_path, bft_path, gatk_path,
              inp_path, ref_path, adp_path, bed_path, out_dir, aligner,
              pic_path, voi_path, java_path, sra_path, verbose, threads, purge):
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
    prep = Prepper(inp_path, out_dir, sra_path).prepInputs()
    samples, sra_list, files = list(), list(), list()
    logger.info('Running MaRS on {0} experiments'.format(len(prep)))
    #summary = Summary(ref_path, bed_path, voi_path, out_dir)
    #samples = config.keys()
    pools = Pool(threads)
    for sample in prep:
        samples.append(prep[sample].sample)
        files.append(prep[sample].files)
        sra_list.append(prep[sample].sra)
    #rone_list = list()
    #rtwo_list = list()
    #name_list = list()
    #for samples in config:
    #    name_list.append(config[samples].sample)
    #    rone_list.append(config[samples].files[0])
    #    rtwo_list.append(config[samples].files[1])

    #sra_list = files
    vcf_list = pools.map(main, zip(repeat(bbduk_path), repeat(aligner_path),
                repeat(smt_path), repeat(bft_path), repeat(gatk_path),
                samples, files, repeat(ref_path), repeat(adp_path),
                repeat(bed_path), repeat(out_dir), repeat(aligner),
                repeat(pic_path), repeat(voi_path),
                repeat(java_path), repeat(sra_path), repeat(purge), sra_list))
    logger.info('Summarizing variant calls from all {0} experiments'.format(len(prep)))
    summary = Summary(ref_path, bed_path, voi_path, out_dir)
    #Sumarize variants of intrest
    summary.getSummary()
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
    ### comment: add path through argument
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
    parser.add_argument('--threads', dest='threads', type=int, default=5,
                        help='Number of threads')
    parser.add_argument('--verbose', action='store_true', 
                        help='Increase verbosity of log file')                        
    parser.add_argument('--purge', action='store_true', 
                        help='Remove intermiediate Fastq and alignment files')                        
    args = parser.parse_args()

    #Validate parsed arguments
    if args.aligner_path is None:
        args.aligner_path = aligner_def[args.aligner]

    if not os.path.exists(args.out_path):
        os.mkdir(args.out_path)

    #single sample experiment.

    #Check if the run command is for batch mode analysis or single sample
    #analysis.
    #If inp_path is empty and rone_path is not, then the experiment is a
    #single sample experiment.
    status = marsBatch(args.bbduk_path, args.aligner_path, args.smt_path,
                args.bft_path, args.gatk_path, args.inp_path, args.ref_path,
                args.adp_path, args.bed_path, args.out_path, args.aligner,
                args.pic_path, args.voi_path, java_def, sra_def, args.verbose, 
                args.threads, args.purge)


####comment NeST overall:

###BWA-MEM: High-quality queries, 
###low-divergent sequences against a large reference genome, 
###Illumina sequence reads to 70bp to 100 bp (1Mbp)

###Bbduk: combine common data quality for trimming, filtering, and mask
### many trimming, filtering, 

###Bowtie: Good for long reference squences
###algining 100s or 1000s of characters (especially for long mammalian genomes)

###Freebayes:Bayesian genetic variant detector
###find small polymorphisms, SNP, indels, MNPs, and complex events
###halotype based(calls variants based on literal sequences of reads)
###detect variants based on alginments


###SNAP:New sequence aligner 3-20x faster
### supports a rich error model 
### good for cheaply match reads
###reads BAM,FASTQ, SAM or BAM

###BBmap:algin reads from all major platforms
### good for highly mutated genomes or reads with long indels
### even covers 100kbp long whole-gene deletions
### metagenome with over 200 million contigs
###large array options, quality control of libraries and sequencing runs
###----BBsplit is also good for binning or refining metagenomic reads