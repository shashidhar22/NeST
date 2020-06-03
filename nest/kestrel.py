import os
import re
import sys
import csv
import glob
import argparse
import subprocess
from multiprocessing import Pool
from itertools import repeat

###reference: https://github.com/paudano/kestrel
###comment: Mapping-free variant caller for short-read Illumina data
###Kestrel is a variant-caller for short-read Illumina data. It does not align 
###sequence reads or perform a de novo assembly. Instead, it breaks reads into k-mers, 
###which are short overlapping fragments of uniform length, and it counts the number of 
###occurences of each unique k-mer in the sequence data. It also finds the ordered k-mers 
###of a reference. With these two k-mer sets, Kestrel searches for patterns of variation. 
###By using a novel local assembly method guided by the reference, it builds one or more 
###sequences over the altered region and call variants from its alignment.
class KestrelVar:

    def __init__(self, rone_path, rtwo_path, ref_path, kanalyze_path, kestrel_path, out_path):
        self.rone_path = rone_path
        self.rtwo_path = rtwo_path
        self.ref_path = os.path.abspath(ref_path)
        self.out_path = os.path.abspath(out_path)
        self.kestrel_path = os.path.abspath(kestrel_path)
        self.kanalyze_path = os.path.abspath(kanalyze_path)
        if not os.path.exists(self.out_path):
            os.mkdir(self.out_path)

        return


    ###comment: running kestrel function
    def run_kestrel(self):

        ikc_path = '{0}/kmertable.ikc'.format(self.out_path)

        ###comment: input command to run the kestrel
        kcmd =['java', '-jar', self.kanalyze_path, 'count', '-k', '31', '-m', 'ikc',
            '--minsize', '15','-o', ikc_path, self.rone_path, self.rtwo_path]

        if not os.path.exists(ikc_path):

            krun = subprocess.Popen(kcmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)
            kout = krun.communicate()[0]
            kret = krun.returncode

            if kret != 0:
                print('Kanalyze failed')
                print(' '.join(kcmd))
                return(ikc_path, None)
#        print(' '.join(kcmd))
        var_path = '{0}/vairants_kes.vcf'.format(self.out_path)
        kcmd = ['java', '-jar', self.kestrel_path, '-r',
            self.ref_path, '-m', 'vcf', '--noanchorboth',
            '--varfilter=coverage:0.5,5', '--loglevel=all',
            '--logfile={}/kestrel.log'.format(self.out_path), '-o', var_path,
            ikc_path]

        krun = subprocess.Popen(kcmd, stdout=subprocess.PIPE,
                                            stderr=subprocess.PIPE, shell=False)
        kout = krun.communicate()[0]
        kret = krun.returncode

#        print(' '.join(kcmd))
        if kret != 0:
            print('Kestrel failed')
            print(' '.join(kcmd))
        return(var_path, kret)


###comment: function for running the kesterl class and the run_kestrel function inside the kestrel class
def kes_runner(rone_path, rtwo_path, ref_path, samtools_path, kanalyze_path, kestrel_path, out_path):
    rone_path = os.path.abspath(rone_path)
    rtwo_path = os.path.abspath(rtwo_path)
    ref_path = os.path.abspath(ref_path)
    kanalyze_path = os.path.abspath(kanalyze_path)
    kestrel_path = os.path.abspath(kestrel_path)
    out_path = os.path.abspath(out_path)
    if not os.path.exists(out_path):
        os.mkdir(out_path)

    kes_run = KestrelVar(rone_path, rtwo_path, ref_path, samtools_path, kanalyze_path, kestrel_path, out_path)
    var_path, kret = kes_run.run_kestrel()
    return(kret, var_path)


###arguments necessary to run the kestrel
if __name__ == '__main__':

     parser = argparse.ArgumentParser(prog='Kestrel E.coli test')
     parser.add_argument('--input', type=str, help="Path to directory with assemblies")
     parser.add_argument('--ref', type=str, help="Path to reference fasta file")
     parser.add_argument('--bwa', type=str, help="Path to bwa aligner")
     parser.add_argument('--sam', type=str, help="Path to samtools")
     parser.add_argument('--kan', type=str, help="Path to kanlyze")
     parser.add_argument('--kes', type=str, help="Path to kestrel")
     parser.add_argument('--out', type=str, help="Path to output directory")

     args = parser.parse_args()
     runner(args.input, args.ref, args.bwa, args.sam, args.kan, args.kes, args.out)
