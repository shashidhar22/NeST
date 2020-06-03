import os
import sys
import time
import logging
import argparse
import subprocess


###reference: https://www.webopedia.com/TERM/O/object_oriented_programming_OOP.html
###Encapsulation: The process of combining elements to create a new entity. A procedure is a type of encapsulation because it combines a series of computer instructions.
###Information hiding: The process of hiding details of an object or function. Information hiding is a powerful programming technique because it reduces complexity.
###Inheritance: a feature that represents the "is a" relationship between different classes.
###Interface: the languages and codes that the applications use to communicate with each other and with the hardware.
###Messaging: Message passing is a form of communication used in parallel programming and object-oriented programming.
###Object: a self-contained entity that consists of both data and procedures to manipulate the data.
###Polymorphism: A programming language's ability to process objects differently depending on their data type or class.
###Procedure: a section of a program that performs a specific task.
###comment: Class: A category of objects. The class defines all the common 
####properties of the different objects that belong to it.
####This python file contains aligners but they have different object/properties
###comment: class BWA alginer 

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
        ###set input paths
        self.bwa_path = bwa_path
        self.out_path = '{0}/alignments'.format(out_path)
        self.ref_path = ref_path
        self.logger = logging.getLogger('NeST.alignment')
        ###if paths don't exist make directories for those paths
        if not os.path.exists(self.out_path):
            os.mkdir(self.out_path)
        #Build index if its not present
        if not os.path.exists('{0}.sa'.format(self.ref_path)):
            self.logger.debug(('Reference fasta file not indexed;'
                'Creating BWA index files'))
            icmd = ['bwa', 'index', ref_path]
            irun = subprocess.Popen(icmd, shell=False, stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL)
            irun.wait()
            if irun.returncode != 0:
                self.logger.error('Reference fasta could not be indexed')
                self.logger.error(' '.join(icmd))
            else:
                self.logger.debug('Reference fasta indexed')
        return

    ###reference: http://bio-bwa.sourceforge.net/bwa.shtml
    ###comment: BWA-MEM, which is the latest, is generally recommended for high-quality queries as it is faster and more accurate. 
    ###BWA-MEM also has better performance than BWA-backtrack for 70-100bp Illumina reads. 
    ###The BWA-MEM algorithm performs local alignment. It may produce multiple primary alignments 
    ###for different part of a query sequence. This is a crucial feature for long sequences. However, 
    ###some tools such as Picard’s markDuplicates does not work with split alignments. One may consider 
    ###to use option -M to flag shorter split hits as secondary.
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
        ###comment: input commands needed to run the bwa
        bwcmd = [self.bwa_path, 'mem', '-t', '4', self.ref_path,
                rone_path, rtwo_path, '>', sam_path]
        ###comment: running bwa with given inputs from above
        bwrun = subprocess.Popen(' '.join(bwcmd), stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL, shell=True)
        bwrun.wait()
        return(sam_path, bwrun.returncode)

###reference: http://bowtie-bio.sourceforge.net/manual.shtml
###comment: Bowtie is an ultrafast, memory-efficient short read aligner geared 
###toward quickly aligning large sets of short DNA sequences (reads) to large genomes. 
###It aligns 35-base-pair reads to the human genome at a rate of 25 million reads per hour on a 
###typical workstation. Bowtie indexes the genome with a Burrows-Wheeler index to keep its memory 
###footprint small: for the human genome, the index is typically about 2.2 GB (for unpaired alignment) 
###or 2.9 GB (for paired-end or colorspace alignment). Multiple processors can be used simultaneously 
###to achieve greater alignment speed. Bowtie can also output alignments in the standard SAM format, 
###allowing Bowtie to interoperate with other tools supporting SAM, including the SAMtools consensus, 
###SNP, and indel callers. Bowtie runs on the command line under Windows, Mac OS X, Linux, and Solaris.
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
        self.logger = logging.getLogger('NeST.alingment')
        if not os.path.exists(self.out_path):
            os.mkdir(self.out_path)
        #Build index if its not present
        if not os.path.exists('{0}.1.bt2'.format(self.ref_path)):
            self.logger.debug(('Reference fasta file not indexed;'
                'Creating Bowtie2 index files'))
            icmd = ['bowtie2-build', ref_path, ref_path]
            irun = subprocess.Popen(icmd, shell=False, stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL)
            irun.wait()
            if irun.returncode != 0:
                self.logger.error('Reference fasta could not be indexed')
                self.logger.error(' '.join(icmd))
            else:
                self.logger.debug('Reference fasta indexed')

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
        ###command to run the bowtie, the inputs needed
        bwcmd = [self.bowtie_path, '-x', self.ref_path, '-1',
                rone_path, '-2', rtwo_path, '--very-sensitive',
                '-p', '4', '-S', sam_path]
        ###command to run the bowtie with given inputs from above
        bwrun = subprocess.Popen(bwcmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=False)
        bwrun.wait()
        if bwrun.returncode:
            self.logger.error(' '.join(bwcmd))

        return(sam_path, bwrun.returncode)


###reference: https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/
###comment: BBMap class
###comment: BBMap is a splice-aware global aligner for DNA and RNA sequencing reads. 
###It can align reads from all major platforms – Illumina, 454, Sanger, Ion Torrent, 
###Pac Bio, and Nanopore. BBMap is fast and extremely accurate, particularly with highly 
###mutated genomes or reads with long indels, even whole-gene deletions over 100kbp long. 
###It has no upper limit to genome size or number of contigs, and has been successfully used 
###for mapping to an 85 gigabase soil metagenome with over 200 million contigs. Additionally, 
###the indexing phase is very fast compared to other aligners.
###BBMap has a large array of options, described in its shell script. It can output 
###many different statistics files, such as an empirical read quality histogram, 
###insert-size distribution, and genome coverage, with or without generating a sam file. 
###As a result, it is useful in quality control of libraries and sequencing runs, or 
###evaluating new sequencing platforms. The derivative program BBSplit is also useful in 
###binning or refining metagenomic reads.
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
        self.logger = logging.getLogger('NeST.alingment')
        if not os.path.exists(self.out_path):
            os.mkdir(self.out_path)
        #Build index if its not present
        if not os.path.exists('{0}.sa'.format(self.out_path)):
            self.logger.debug(('Reference fasta file not indexed;'
                'Creating BBMap index files'))
            icmd = ['bbmap.sh', 'ref={0}'.format(ref_path)]
            irun = subprocess.Popen(icmd, shell=False, stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL)
            irun.wait()
            if irun.returncode != 0:
                self.logger.error('Reference fasta could not be indexed')
            else:
                self.logger.debug('Reference fasta indexed')

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
        ###comment: input commands for bbmap
        bbcmd = [self.bbmap_path, 'ref={0}'.format(self.ref_path),
                'in={0}'.format(rone_path), 'in2={0}'.format(rtwo_path),
                'out={0}'.format(sam_path)]
        bbrun = subprocess.Popen(bbcmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=False)
        bbrun.wait()

        return(sam_path, bbrun.returncode)

###refernece: https://bioinformaticshome.com/tools/descriptions/SNAP.html
###comment: A gene prediction tool loosely based on HMM. It has functionality to generate graphical plots.
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
        self.logger = logging.getLogger('NeST.alingment')
        if not os.path.exists(self.out_path):
            os.mkdir(self.out_path)
        #Build index if its not present
        if not os.path.exists('{0}.sa'.format(self.out_path)):
            ref_out = os.path.dirname(self.ref_path)
            self.logger.debug(('Reference fasta file not indexed;'
                'Creating SNAP index files'))
            icmd = ['snap-aligner', 'index', ref_path, ref_out]
            irun = subprocess.Popen(icmd, shell=False, stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL)
            irun.wait()
            if irun.returncode != 0:
                self.logger.error('Reference fasta could not be indexed')
            else:
                self.logger.debug('Reference fasta indexed')


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
