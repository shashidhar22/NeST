import os
import sys
import time
import logging
import argparse
import subprocess


###reference: https://gatk.broadinstitute.org/hc/en-us
###comment: The GATK is the industry standard for identifying SNPs and indels in germline DNA and RNAseq 
###data. Its scope is now expanding to include somatic short variant calling, and to tackle copy number 
###(CNV) and structural variation (SV). In addition to the variant callers themselves, the GATK also 
###includes many utilities to perform related tasks such as processing and quality control of 
###high-throughput sequencing data, and bundles the popular Picard toolkit.
###comment: These tools were primarily designed to process exomes and whole genomes generated with 
###Illumina sequencing technology, but they can be adapted to handle a variety of other technologies 
###and experimental designs. And although it was originally developed for human genetics, the GATK has 
###since evolved to handle genome data from any organism, with any level of ploidy.
class GenAnTK:

    def __init__(self, gatk_path, out_path, java, picard_path):
        self.gatk_path = gatk_path
        self.out_path = out_path
        self.java = java
        self.picard_path = picard_path

        return
    ###reference: https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller
    ###comment: The HaplotypeCaller is capable of calling SNPs and indels simultaneously 
    ###via local de-novo assembly of haplotypes in an active region. In other words, 
    ###whenever the program encounters a region showing signs of variation, it discards the 
    ###existing mapping information and completely reassembles the reads in that region. This 
    ###allows the HaplotypeCaller to be more accurate when calling regions that are traditionally 
    ###difficult to call, for example when they contain different types of variants close to each 
    ###other. It also makes the HaplotypeCaller much better at calling indels than position-based 
    ###callers like UnifiedGenotyper.
    def hapCaller(self, bam_path, ref_path, sam_name):
        vcf_path = '{0}/{1}_variants_gatk.vcf'.format(self.out_path, sam_name)
        ref_base = os.path.splitext(ref_path)[0]
        dict_file = '{0}.dict'.format(ref_base)
        if not os.path.exists(dict_file):
            pcmd = [self.picard_path, 'CreateSequenceDictionary', 'R={0}'.format(ref_path), 
                    'O={0}'.format(dict_file)]
            prun = subprocess.Popen(pcmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=False)
            prun.wait()
            
        hcmd = [self.gatk_path, 'HaplotypeCaller', '-R', ref_path,
            '-I', bam_path, '-O', vcf_path, '--sample-name', sam_name,
            '--genotyping-mode', 'DISCOVERY', '--output-mode', 'EMIT_VARIANTS_ONLY']
        hrun = subprocess.Popen(hcmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=False)
        hrun.wait()
        if hrun.returncode != 0:
            print(' '.join(hcmd))
        return(vcf_path, hrun.returncode)

###reference: https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionVI/lessons/02_variant-calling.html
###comment: FreeBayes is a Bayesian genetic variant detector designed to find small polymorphisms, 
###specifically SNPs (single-nucleotide polymorphisms), indels (insertions and deletions), MNPs 
###(multi-nucleotide polymorphisms), and complex events (composite insertion and substitution events) 
###smaller than the length of a short-read sequencing alignment
###FreeBayes is haplotype-based, in the sense that it calls variants based on the literal sequences of reads 
###aligned to a particular target, not their precise alignment. This model is a straightforward 
###generalization of previous ones (e.g. PolyBayes, samtools, GATK) which detect or report variants 
###based on alignments. This method avoids one of the core problems with alignment-based variant 
###detection— that identical sequences may have multiple possible alignments:
class FreeBayes:

    def __init__(self, free_path, out_path):
        self.free_path = free_path
        self.out_path = out_path
        return
    def freeBayes(self, bam_path, ref_path, sam_name):
        vcf_path = '{0}/{1}_variants_freebayes.vcf'.format(self.out_path, sam_name)
        fcmd = '{0} -f {1} {2} > {3}'.format(self.free_path, ref_path, bam_path, vcf_path)
        frun = subprocess.Popen(fcmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True)
        frun.wait()
        if frun.returncode != 0:
            print(' '.join(fcmd))
        return(vcf_path, frun.returncode)

###reference: https://broadinstitute.github.io/picard/
###comment: Picard is a set of command line tools for manipulating high-throughput sequencing (HTS) 
###data and formats such as SAM/BAM/CRAM and VCF. These file formats are defined in the Hts-specs 
###repository. See especially the SAM specification and the VCF specification.
class Picard:

    def __init__(self, java, pic_path, out_path):
        self.java = java
        self.pic_path = pic_path
        self.out_path = out_path
        self.logger = logging.getLogger('NeST.picard')

    def picard(self, bam_path, sam_name):
        add_path = '{0}/alignments/{1}_RG.bam'.format(self.out_path, os.path.splitext(os.path.basename(bam_path))[0])
        acmd = [self.pic_path, 'AddOrReplaceReadGroups',
            'I='+bam_path , 'O='+add_path, 'SORT_ORDER=coordinate', 'RGID={0}'.format(sam_name),
            'RGLB=ExomeSeq', 'RGPL=Illumina', 'RGPU=HiSeq2500', 'RGSM={0}'.format(sam_name),
            'RGCN=AtlantaGenomeCenter', 'RGDS=ExomeSeq', 'RGDT=2016-08-24', 'RGPI=null',
            'RGPG={0}'.format(sam_name), 'RGPM={0}'.format(sam_name), 'CREATE_INDEX=true']
        arun = subprocess.Popen(acmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=False)
        arun.wait()

        if arun.returncode != 0:
            self.logger.error('Picard add read group information failed running the following command : {0}'.format(' '.join(acmd)))
        return(add_path, arun.returncode)
