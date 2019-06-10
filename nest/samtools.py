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
        self.logger = logging.getLogger('NeST.samtools')

    def fixmate(self, bam_path):
        base = os.path.splitext(os.path.basename(bam_path))[0]
        obam_path = '{0}/alignments/{1}_FM.bam'.format(self.out_path, base)
        fmcmd = [self.sam_path, 'fixmate', '-O', 'BAM',
                bam_path, obam_path]
        fmrun = subprocess.Popen(fmcmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=False)
        fmrun.wait()
        if fmrun.returncode != 0:
            self.logger.error('Samtools fixmate failed running the following command : {0}'.format(' '.join(fmcmd)))
            print(' '.format(fmcmd))

        return(obam_path, fmrun.returncode)


    def sort(self, bam_path):
        base = os.path.splitext(os.path.basename(bam_path))[0]
        obam_path = '{0}/alignments/{1}_SR.bam'.format(self.out_path, base)
        stcmd = [self.sam_path, 'sort', '-O', 'BAM',
                '-o', obam_path, bam_path]
        strun = subprocess.Popen(stcmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=False)
        strun.wait()
        if strun.returncode != 0:
            self.logger.error('Samtools sort failed running the following command : {0}'.format(' '.join(stcmd)))

        return(obam_path, strun.returncode)

    def dedup(self, bam_path):
        base = os.path.splitext(os.path.basename(bam_path))[0]
        obam_path = '{0}/alignments/{1}_DD.bam'.format(self.out_path, base)
        ddcmd = [self.sam_path, 'rmdup', bam_path, obam_path]
        ddrun = subprocess.Popen(ddcmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=False)
        ddrun.wait()
        if ddrun.returncode != 0:
            self.logger.error('Samtools sort failed running the following command : {0}'.format(' '.join(ddcmd)))

        return(obam_path, ddrun.returncode)

    def addreadgroup(self, bam_path, sam_name):
        run_time = time.strftime('%d/%m/%Y')
        base = os.path.splitext(os.path.basename(bam_path))[0]
        abam_path = '{0}/alignments/{1}_RG.bam'.format(self.out_path, base)
        acmd = [self.sam_path, 'addreplacerg', '-O', 'BAM',
                '-r', '"@RG\tID:{0}\tPG:{0}\tSM:{0}\tPM:{0}\tLB=MaRS\tPL=Illumina\tPU=MiSeq\tDT={1}\tPI=null"'.format(sam_name, run_time),
                '-o', abam_path, bam_path]
        arun = subprocess.Popen(acmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=False)
        arun.wait()
        icmd = [self.sam_path, 'index', abam_path]
        irun = subprocess.Popen(icmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=False)
        irun.wait()
        if arun.returncode != 0:
            print(' '.join(acmd))
        return(abam_path, arun.returncode)

    def pileup(self, ref_path, bam_path, sam_name):
        obcf_path = '{0}/{1}_variants.bcf'.format(self.out_path, sam_name)
        fai_path = '{0}.fai'.format(ref_path)
        if not os.path.exists(fai_path):
            facmd = [self.sam_path, 'faidx', ref_path]
            farun = subprocess.Popen(facmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=False)
        mpcmd = [self.sam_path, 'mpileup', '-go', obcf_path,
                '-f', ref_path, bam_path]
        mprun = subprocess.Popen(mpcmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=False)
        mprun.wait()
        if mprun.returncode != 0:
            self.logger.error('Samtools mpileup failed running the following command : {0}'.format(' '.join(mpcmd)))

        return(obcf_path, mprun.returncode)

    def bcfindex(self, bcf_path):

        bicmd = [self.bft_path, 'index', bcf_path]
        birun = subprocess.Popen(bicmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=False)
        birun.wait()
        if birun.returncode != 0:
            self.logger.error('Bcftools index failed running the following command : {0}'.format(' '.join(bicmd)))

        return(birun.returncode)

    def bcftools(self, bcf_path, sam_name):
        ovcf_path = '{0}/{1}_variants_samtools.vcf'.format(self.out_path, sam_name)
        btcmd = [self.bft_path, 'call', 
                '--multiallelic-caller', '--variants-only', '-O', 'v',
                '-s', sam_name, '-o', ovcf_path, bcf_path]
        btrun = subprocess.Popen(btcmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=False)
        btrun.wait()
        if btrun.returncode != 0:
            self.logger.error('Bcftools call failed running the following command : {0}'.format(' '.join(btcmd)))

        return(ovcf_path, btrun.returncode)


    def bcfstats(self, vcf_path, ref_path):
        stats_path = '{0}/variants.stats'.format(self.out_path)
        bscmd = [self.bft_path, 'stats', '-F', ref_path, '-s' , '-',
                vcf_path, '>', stats_path]
        bsrun = subprocess.Popen(' '.join(bscmd), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True)
        bsrun.wait()
#        if bsrun.returncode != 0:
#            logger.error('Bcftools stats failed running the following command : {0}'.format(' '.join(bscmd)))

        return(stats_path, bsrun.returncode)
