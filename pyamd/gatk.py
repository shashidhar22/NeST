import os
import sys
import time
import logging
import argparse
import subprocess

logger = logging.getLogger('GATK')
logger.setLevel(logging.ERROR)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s:%(name)s:%(levelname)s:%(message)s')
ch.setFormatter(formatter)
logger.addHandler(ch)


class GenAnTK:

    def __init__(self, gatk_path, out_path, java):
        self.gatk_path = gatk_path
        self.out_path = out_path
        self.java = java

        return

    def hapCaller(self, bam_path, ref_path, sam_name):
        vcf_path = '{0}/{1}_variants_gatk.vcf'.format(self.out_path, sam_name)

        hcmd = [self.java, '-jar', self.gatk_path, '-T', 'HaplotypeCaller', '-R', ref_path,
            '-I', bam_path, '-o', vcf_path, '-nct', '1', '-sn', sam_name,
            '-gt_mode', 'DISCOVERY']
        hrun = subprocess.Popen(hcmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=False)
        hrun.wait()
        if hrun.returncode != 0:
            logger.error('GATK HaplotypeCaller failed running the following command : {0}'.format(' '.join(hcmd)))
        return(vcf_path, hrun.returncode)


class Picard:

    def __init__(self, java, pic_path, out_path):
        self.java = java
        self.pic_path = pic_path
        self.out_path = out_path


    def picard(self, bam_path, sam_name):
        add_path = '{0}/{1}_RG.bam'.format(self.out_path, os.path.splitext(os.path.basename(bam_path))[0])
        acmd = [self.java, '-Xmx1g', '-jar', self.pic_path, 'AddOrReplaceReadGroups',
            'I='+bam_path , 'O='+add_path, 'SORT_ORDER=coordinate', 'RGID={0}'.format(sam_name),
            'RGLB=ExomeSeq', 'RGPL=Illumina', 'RGPU=HiSeq2500', 'RGSM={0}'.format(sam_name),
            'RGCN=AtlantaGenomeCenter', 'RGDS=ExomeSeq', 'RGDT=2016-08-24', 'RGPI=null',
            'RGPG={0}'.format(sam_name), 'RGPM={0}'.format(sam_name), 'CREATE_INDEX=true']
        arun = subprocess.Popen(acmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=False)
        arun.wait()
        if arun.returncode != 0:
            logger.error('Picard AddOrReplaceReadGroups failed running the following command : {0}'.format(' '.join(acmd)))

        return(add_path, arun.returncode)
