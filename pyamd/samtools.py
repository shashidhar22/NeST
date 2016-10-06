import os
import sys
import time
import logging
import argparse
import subprocess


class Samtools:
    
    def __init__(self, sam_path, bft_path, out_path):
        self.sam_path = sam_path
        self.out_path = out_path
        self.bft_path = bft_path
    
    def fixmate(self, bam_path):
        base = os.path.splitext(os.path.basename(bam_path))[0]
        obam_path = '{0}/{1}_fixmate.bam'.format(self.out_path, base)
        fmcmd = [self.sam_path, 'fixmate', '-O', 'BAM',
                bam_path, obam_path]
        print(' '.join(fmcmd))
        fmrun = subprocess.Popen(fmcmd, shell=False)
        fmrun.wait()

        return(obam_path, fmrun.returncode)


    def sort(self, bam_path):
        base = os.path.splitext(os.path.basename(bam_path))[0]
        obam_path = '{0}/{1}_sorted.bam'.format(self.out_path, base)
        stcmd = [self.sam_path, 'sort', '-O', 'bam',
                '-o', obam_path, bam_path]

        strun = subprocess.Popen(stcmd, shell=False)
        strun.wait()

        return(obam_path, strun.returncode)

    def pileup(self, ref_path, bam_path):
        obcf_path = '{0}/variants.bcf'.format(self.out_path)
        
        mpcmd = [self.sam_path, 'mpileup', '-go', obcf_path,
                '-f', ref_path, bam_path]
        
        mprun = subprocess.Popen(mpcmd, shell=False)
        mprun.wait()
        
        return(obcf_path, mprun.returncode)

    def bcfindex(self, bcf_path):
        
        bicmd = [self.bft_path, 'index', bcf_path]
        
        birun = subprocess.Popen(bicmd, shell=False)
        birun.wait()
        
        return(birun.returncode)

    def bcftools(self, bcf_path, bed_path):
        ovcf_path = '{0}/variants.vcf'.format(self.out_path)
        btcmd = [self.bft_path, 'call', '--skip-variants', 'indels',
                '--multiallelic-caller', '--variants-only', '-O', 'v', 
                '-o', ovcf_path, bcf_path]        
#        btcmd = [self.bft_path, 'call', '-V', 'indels', '-vmO', 'v', 
#                '-o', ovcf_path, bcf_path]
#        print(' '.join(btcmd))
        btrun = subprocess.Popen(btcmd, shell=False)
        btrun.wait()
        
        return(ovcf_path, btrun.returncode)
      

    def bcfstats(self, vcf_path, ref_path):
        stats_path = '{0}/variants.stats'.format(self.out_path)
        bscmd = [self.bft_path, 'stats', '-F', ref_path, '-s' , '-',
                vcf_path, '>', stats_path]
        
        bsrun = subprocess.Popen(' '.join(bscmd), shell=True)
        bsrun.wait()

        return(stats_path, bsrun.returncode) 
