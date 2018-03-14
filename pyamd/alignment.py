import os
import sys
import time
import logging
import argparse
import subprocess


class Bwa:
    '''Bwa class runs BWA mem in standard settings to align the samples reads
    against the reference genome
    Attributes:
        1. bwa_path : Path to BWA executable
        2. out_path : Output path for the sample
        3. ref_path : Path to reference Fasta file
    Class variables:
        1. self.bwa_path : Path to BWA mem
        2. self.out_path : Output path for the sample
        3. self.ref_path : Path to reference Fasta file
    '''
    def __init__(self, bwa_path, out_path, ref_path):
        self.bwa_path = bwa_path
        self.out_path = '{0}/alignments'.format(out_path)
        self.ref_path = ref_path
        if not os.path.exists(self.out_path):
            os.mkdir(self.out_path)
        return

    def bwamem(self, rone_path, rtwo_path):
        '''Bwamem methods runs BWA mem using default parameters to align R1 and
        R2 reads to the reference genome.
        Attributes :
            1. rone_path : Path to R1 Fastq file
            2. rtwo_path : Path to R2 Fastq file
        Return values are:
            1. sam_path : Path to sam files
            2. bwrun.returncode : Returncode from BWA execution
        '''
        sam_path = '{0}/output.sam'.format(self.out_path)
        bwcmd = [self.bwa_path, 'mem', '-t', '4', self.ref_path,
                rone_path, rtwo_path, '>', sam_path]
        bwrun = subprocess.Popen(' '.join(bwcmd), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True)
        bwrun.wait()

        return(sam_path, bwrun.returncode)


class Bowtie:
    '''Bowtie class runs Bowtie2 in standard settings to align the samples reads
    against the reference genome
    Attributes:
        1. bowtie_path : Path to Bowtie2 executable
        2. out_path : Output path for the sample
        3. ref_path : Path to reference Fasta file
    Class variables:
        1. self.bowtie_path : Path to Bowtie2
        2. self.out_path : Output path for the sample
        3. self.ref_path : Path to reference Fasta file
    '''
    def __init__(self, bowtie_path, out_path, ref_path):
        self.bowtie_path =  bowtie_path
        self.out_path = '{0}/alignments'.format(out_path)
        self.ref_path = ref_path
        if not os.path.exists(self.out_path):
            os.mkdir(self.out_path)
        return

    def bowtie(self, rone_path, rtwo_path):
        '''Bowtie methods runs Bowtie2 using default parameters to align R1 and
        R2 reads to the reference genome.
        Attributes :
            1. rone_path : Path to R1 Fastq file
            2. rtwo_path : Path to R2 Fastq file
        Return values are:
            1. sam_path : Path to sam files
            2. bwrun.returncode : Returncode from Bowtie2 execution
        '''
        sam_path = '{0}/output.sam'.format(self.out_path)
        bwcmd = [self.bowtie_path, '-x', self.ref_path, '-1',
                rone_path, '-2', rtwo_path, '--very-sensitive',
                '-p', '4', '-S', sam_path]
        bwrun = subprocess.Popen(bwcmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=False)
        bwrun.wait()

        return(sam_path, bwrun.returncode)


class BBMap:
    '''BBMap class runs BBMap in standard settings to align the samples reads
    against the reference genome
    Attributes:
        1. bbmap_path : Path to BBMap executable
        2. out_path : Output path for the sample
        3. ref_path : Path to reference Fasta file
    Class variables:
        1. self.bbmap_path : Path to BBMap
        2. self.out_path : Output path for the sample
        3. self.ref_path : Path to reference Fasta file
    '''
    def __init__(self, bbmap_path, out_path, ref_path):
        self.bbmap_path = bbmap_path
        self.out_path = out_path
        self.ref_path =  ref_path
        if not os.path.exists(self.out_path):
            os.mkdir(self.out_path)

    def bbmap(self, rone_path, rtwo_path):
        '''Bbmap method runs BBMap using default parameters to align R1 and
        R2 reads to the reference genome.
        Attributes :
            1. rone_path : Path to R1 Fastq file
            2. rtwo_path : Path to R2 Fastq file
        Return values are:
            1. sam_path : Path to sam files
            2. bbrun.returncode : Returncode from BBMap execution
        '''
        sam_path = '{0}/output.sam'.format(self.out_path)
        bbcmd = [self.bbmap_path, 'ref={0}'.format(self.ref_path),
                'in={0}'.format(rone_path), 'in2={0}'.format(rtwo_path),
                'out={0}'.format(sam_path)]
        bbrun = subprocess.Popen(bbcmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=False)
        bbrun.wait()

        return(sam_path, bbrun.returncode)

class Snap:
    '''Snap class runs SNAP in standard settings to align the samples reads
    against the reference genome
    Attributes:
        1. snap_path : Path to SNAP executable
        2. out_path : Output path for the sample
        3. ref_path : Path to reference Fasta file
    Class variables:
        1. self.snap_path : Path to SNAP
        2. self.out_path : Output path for the sample
        3. self.ref_path : Path to reference Fasta file
    '''
    def __init__(self, snap_path, out_path, ref_path):
        self.snap_path = snap_path
        self.out_path = out_path
        self.ref_path = os.path.dirname(ref_path)
        if not os.path.exists(self.out_path):
            os.mkdir(self.out_path)


    def snap(self, rone_path, rtwo_path):
        '''Snap method runs SNAP using default parameters to align R1 and
        R2 reads to the reference genome.
        Attributes :
            1. rone_path : Path to R1 Fastq file
            2. rtwo_path : Path to R2 Fastq file
        Return values are:
            1. sam_path : Path to sam files
            2. srun.returncode : Returncode from SNAP execution
        '''
        sam_path = '{0}/output.sam'.format(self.out_path)
        scmd = [self.snap_path, 'paired', self.ref_path, rone_path, rtwo_path,
                '-t', '4', '-o', '-sam', sam_path]
        srun = subprocess.Popen(scmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=False)
        srun.wait()
        return(sam_path, srun.returncode)
