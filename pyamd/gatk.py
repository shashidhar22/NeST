import os
import sys
import time
import logging
import argparse
import subprocess


class GenAnTK:
    
    def __init__(self, gatk_path, pic_path, out_path):
        self.gatk_path = gatk_path
        self.out_path = out_path
        self.pic_path = pic_path
        return

    def picard(self, bam_path):
        add_path = '{0}/{1}_RG.bam'.format(self.out_path, os.path.splitext(os.path.basename(bam_path))[0])
        acmd = ['java', '-Xmx1g', '-jar', self.pic_path, 'AddOrReplaceReadGroups',
            'I='+bam_path , 'O='+add_path, 'SORT_ORDER=coordinate', 'RGID=Test', 
            'RGLB=ExomeSeq', 'RGPL=Illumina', 'RGPU=HiSeq2500', 'RGSM=Test', 
            'RGCN=AtlantaGenomeCenter', 'RGDS=ExomeSeq', 'RGDT=2016-08-24', 'RGPI=null', 
            'RGPG=Test', 'RGPM=Test', 'CREATE_INDEX=true']
    
        arun = subprocess.Popen(acmd, shell=False)
        arun.wait()
        return(add_path, arun.returncode)

    def hapCaller(self, bam_path, ref_path):
        vcf_path = '{0}/variants_gatk.vcf'.format(self.out_path)
 
        hcmd = ['java', '-jar', self.gatk_path, '-T', 'HaplotypeCaller', '-R', ref_path,
            '-I', bam_path, '-o', vcf_path, '-nct', '1', 
            '-gt_mode', 'DISCOVERY']
        print(hcmd)
        hrun = subprocess.Popen(hcmd, shell=False)
        hrun.wait()
        return(vcf_path, hrun.returncode)
