import os
import sys
import time
import logging
import argparse
import subprocess



class GenAnTK:

    def __init__(self, gatk_path, out_path, java, picard_path):
        self.gatk_path = gatk_path
        self.out_path = out_path
        self.java = java
        self.picard_path = picard_path

        return

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
