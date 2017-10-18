import os
import sys
import time
import logging
import argparse
import subprocess

logger = logging.getLogger('Alignment')
logger.setLevel(logging.ERROR)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s:%(name)s:%(levelname)s:%(message)s')
ch.setFormatter(formatter)
logger.addHandler(ch)

class Bwa:

    def __init__(self, bwa_path, out_path, ref_path):
        self.bwa_path = bwa_path
        self.out_path = '{0}/alignments'.format(out_path)
        self.ref_path = ref_path
        self.logger = logging.getLogger('Mars.sample_runner.BWA')
        if not os.path.exists(self.out_path):
            os.mkdir(self.out_path)
        return

    def bwamem(self, rone_path, rtwo_path):
        sam_path = '{0}/output.sam'.format(self.out_path)
        bwcmd = [self.bwa_path, 'mem', '-t', '4', self.ref_path,
                rone_path, rtwo_path, '>', sam_path]
        bwrun = subprocess.Popen(' '.join(bwcmd), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True)
        bwrun.wait()
        if bwrun.returncode != 0:
            self.logger.error('BWA failed running the following command : {0}'.format(' '.join(bwcmd)))

        return(sam_path, bwrun.returncode)


class Bowtie:

    def __init__(self, bowtie_path, out_path, ref_path):
        self.bowtie_path =  bowtie_path
        self.out_path = '{0}/alignments'.format(out_path)
        self.ref_path = ref_path
        self.logger = logging.getLogger('Mars.sample_runner.Bowtie2')
        if not os.path.exists(self.out_path):
            os.mkdir(self.out_path)
        return

    def bowtie(self, rone_path, rtwo_path):
        sam_path = '{0}/output.sam'.format(self.out_path)
        bwcmd = [self.bowtie_path, '-x', self.ref_path, '-1',
                rone_path, '-2', rtwo_path, '--very-sensitive',
                '-p', '4', '-S', sam_path]
        bwrun = subprocess.Popen(bwcmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=False)
        bwrun.wait()
        if bwrun.returncode != 0:
            self.logger.error('Bowtie2 failed running the following command : {0}'.format(' '.join(bwcmd)))

        return(sam_path, bwrun.returncode)


class BBMap:
    def __init__(self, bbmap_path, out_path, ref_path):
        self.bbmap_path = bbmap_path
        self.out_path = out_path
        self.ref_path =  ref_path
        self.logger = logging.getLogger('Mars.sample_runner.BBMap')
        if not os.path.exists(self.out_path):
            os.mkdir(self.out_path)

    def bbmap(self, rone_path, rtwo_path):
        sam_path = '{0}/output.sam'.format(self.out_path)
        bbcmd = [self.bbmap_path, 'ref={0}'.format(self.ref_path),
                'in={0}'.format(rone_path), 'in2={0}'.format(rtwo_path),
                'out={0}'.format(sam_path)]
        bbrun = subprocess.Popen(bbcmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=False)
        bbrun.wait()
        if bbrun.returncode != 0:
            self.logger.error('BBMap failed running the following command : {0}'.format(' '.join(bbcmd)))

        return(sam_path, bbrun.returncode)

class Snap:
    def __init__(self, snap_path, out_path, ref_path):
        self.snap_path = snap_path
        self.out_path = out_path
        self.ref_path = os.path.dirname(ref_path)
        self.logger = logging.getLogger('Mars.sample_runner.Snap')
        if not os.path.exists(self.out_path):
            os.mkdir(self.out_path)


    def snap(self, rone_path, rtwo_path):
        sam_path = '{0}/output.sam'.format(self.out_path)
        scmd = [self.snap_path, 'paired', self.ref_path, rone_path, rtwo_path,
                '-t', '4', '-o', '-sam', sam_path]
        srun = subprocess.Popen(scmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=False)
        srun.wait()
        if srun.returncode != 0:
            self.logger.error('Snap failed running the following command : {0}'.format(' '.join(scmd)))

        return(sam_path, srun.returncode)
