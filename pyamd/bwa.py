import os
import sys
import time
import logging
import argparse
import subprocess


class Bwa:
    
    def __init__(self, bwa_path, out_path, ref_path):
        self.bwa_path = bwa_path
        self.out_path = '{0}/alignments'.format(out_path)
        self.ref_path = ref_path
        
        if not os.path.exists(self.out_path):
            os.mkdir(self.out_path)
        return

    def bwamem(self, rone_path, rtwo_path):
        sam_path = '{0}/output.sam'.format(self.out_path)
        bwcmd = [self.bwa_path, 'mem', '-t', '4', self.ref_path,
                rone_path, rtwo_path, '>', sam_path]
        print(' '.join(bwcmd))
        bwrun = subprocess.Popen(' '.join(bwcmd), shell=True)
        bwrun.wait()

        return(sam_path, bwrun.returncode)



