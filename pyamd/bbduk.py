import os
import sys
import time
import logging
import argparse
import subprocess


class QualCheck:

    def __init__(self, bbduk_path, adp_path, out_path):
        self.bbduk_path = bbduk_path
        self.adp_path = adp_path
        self.out_path = '{0}/CleanedFastq'.format(out_path)
        self.logger = logging.getLogger('Mars.sample_runner.BBDuk')
        if not os.path.exists(self.out_path):
            os.mkdir(self.out_path)


    def bbduk(self, rone_path, rtwo_path):
#java -ea -Xmx6000m -cp .../current jgi.BBDukF ktrim=r k=27 hdist=1 edist=0 mink=4 ref=adapters.fa qtrim=rl trimq=30 minlength=50 qin=33 in=input1.fastq in2=input2.fastq out=output1.fastq out2=output2.fastq
        #Setup output paths
        brone = os.path.splitext(os.path.basename(rone_path))[0]
        brtwo = os.path.splitext(os.path.basename(rtwo_path))[0]
        orone_path = '{0}/{1}_cleaned.fq'.format(self.out_path, brone)
        ortwo_path = '{0}/{1}_cleaned.fq'.format(self.out_path, brtwo)
        stats_path = '{0}/{1}_stats.txt'.format(self.out_path, brone)
        #Set up the command
        bbcmd = [self.bbduk_path, '-Xmx1g', 'k=27', 'hdist=1', 'edist=0', 'ktrim=l',
                'mink=4', 'ref={0}'.format(self.adp_path), 'qtrim=rl',
                'trimq=30', 'minlength=50', 'qin=33', 'overwrite=true',
                'in={0}'.format(rone_path), 'in2={0}'.format(rtwo_path),
                'out={0}'.format(orone_path), 'out2={0}'.format(ortwo_path),
                'stats={0}'.format(stats_path)]

        #Run bbduk
        bbrun = subprocess.Popen(bbcmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=False)
        bbrun.wait()
        if bbrun.returncode != 0:
            self.logger.error('BBDuk failed running the following command : {0}'.format(' '.join(bbcmd)))
        return(orone_path, ortwo_path, bbrun.returncode)
