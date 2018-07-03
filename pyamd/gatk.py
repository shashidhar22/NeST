import os
import sys
import time
import logging
import argparse
import subprocess



class GenAnTK:

    def __init__(self, gatk_path, out_path, java):
        self.gatk_path = gatk_path
        self.out_path = out_path
        self.java = java

        return

    def hapCaller(self, bam_path, ref_path, sam_name):
        vcf_path = '{0}/{1}_variants_gatk.vcf'.format(self.out_path, sam_name)

        hcmd = [self.gatk_path, '-T', 'HaplotypeCaller', '-R', ref_path,
            '-I', bam_path, '-o', vcf_path, '-nct', '1', '-sn', sam_name,
            '-gt_mode', 'DISCOVERY' ]
        hrun = subprocess.Popen(hcmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=False)
        hrun.wait()
        if hrun.returncode != 0:
            print(' '.join(hcmd))
        return(vcf_path, hrun.returncode)


class Picard:

    def __init__(self, java, pic_path, out_path):
        self.java = java
        self.pic_path = pic_path
        self.out_path = out_path


    def fixmate(self, bam_path, sam_name):
        fix_path = '{0}/{1}_FM_SO.bam'.format(self.out_path, os.path.splitext(os.path.basename(bam_path))[0])
        fcmd = [self.pic_path, 'FixMateInformation', 'SO=coordinate', 'AS=false', 'I={0}'.format(bam_path),
                'O={0}'.format(fix_path)]
        frun = subprocess.Popen(fcmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=False)
        frun.wait()

        if frun.returncode != 0:
            print(' '.join(fcmd))

        return(fix_path, frun.returncode)
    

    def picard(self, bam_path, sam_name):
        add_path = '{0}/{1}_RG.bam'.format(self.out_path, os.path.splitext(os.path.basename(bam_path))[0])
        acmd = [self.pic_path, 'AddOrReplaceReadGroups',
            'I='+bam_path , 'O='+add_path, 'SORT_ORDER=coordinate', 'RGID={0}'.format(sam_name),
            'RGLB=ExomeSeq', 'RGPL=Illumina', 'RGPU=HiSeq2500', 'RGSM={0}'.format(sam_name),
            'RGCN=AtlantaGenomeCenter', 'RGDS=ExomeSeq', 'RGDT=2016-08-24', 'RGPI=null',
            'RGPG={0}'.format(sam_name), 'RGPM={0}'.format(sam_name), 'CREATE_INDEX=true',
            'VALIDATION_STRINGENCY=LENIENT']
        arun = subprocess.Popen(acmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=False)
        arun.wait()

        if arun.returncode != 0:
            print(' '.join(acmd))
        return(add_path, arun.returncode)
