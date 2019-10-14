import os
import sys
import time
import logging
import argparse
import subprocess



class QualCheck:
    '''QualCheck class is written to filter reads from a sample, based on
    adapter contamination and low quality reads using BBDuk. The reads are
    trimmed from both sides, reads are scanned across using kmers of max length
    27, and minimum length 4. hdist specfies a hamming distance of 1, that is,
    a max of one mismatch between the kmer and the adapter sequences.
    Sections with average quality less 30 are trimmed. Reads smaller than 50
    are excluded from the study.
    Attributes:
        1. bbduk_path : Path to BBDuk executable
        2. adp_path   : Path to Fasta file with adapter sequences that need to
                        be trimmed
        3. out_path   : Output path for the sample
        4. java : Java executable path
    Class variables:
        1. self.bbduk_path : Path to BBDuk
        2. self.adp_path : Path to Fasta file with adapter sequences that need
                           to be trimmed
        3. self.java : Java executable path
        4. self.out_path : Output path for the sample
    '''
    def __init__(self, bbduk_path, adp_path, out_path, java):
        self.bbduk_path = bbduk_path
        self.adp_path = adp_path
        self.java = java
        self.out_path = '{0}/CleanedFastq'.format(out_path)

        if not os.path.exists(self.out_path):
            os.mkdir(self.out_path)


    def bbduk(self, rone_path, rtwo_path):
        '''BBduk method runs the bbduk on R1 and R2 of the sample. The trimming
        parameters to the following values:
            1. ktrim=r          4. edist=0          7. trimq=30
            2. k=27             5. mink=default     8. minlength=50
            3. hdist=1          6. qtrim=rl         9. qin=33
        The method has the following return values:
            1. orone_path : Filtered R1 file
            2. ortwo_path : Filtered R2 file
            3. bbrun.returncode : BBDuk execution returncode
        '''
        brone = os.path.splitext(os.path.basename(rone_path))[0]
        brtwo = os.path.splitext(os.path.basename(rtwo_path))[0]
        orone_path = '{0}/{1}_cleaned.fq'.format(self.out_path, brone)
        ortwo_path = '{0}/{1}_cleaned.fq'.format(self.out_path, brtwo)
        stats_path = '{0}/{1}_stats.txt'.format(self.out_path, brone)
        #Set up the command
        #Change trimq to 30 after insilico exp
        bbcmd = [self.bbduk_path, '-Xmx1g', 'k=27', 'hdist=1', 'edist=1', 'ktrim=l',
                'ref={0}'.format(self.adp_path), 'qtrim=rl', 'minlength=50',
                'trimq=20', 'qin=33', 'overwrite=true', 'mink=4',
                'in={0}'.format(rone_path), 'in2={0}'.format(rtwo_path),
                'out={0}'.format(orone_path), 'out2={0}'.format(ortwo_path),
                'stats={0}'.format(stats_path)]

        #Run bbduk
        bbrun = subprocess.Popen(bbcmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=False)
        bbrun.wait()
        if bbrun.returncode != 0:
            print(' '.join(bbcmd))
        return(orone_path, ortwo_path, bbrun.returncode)
