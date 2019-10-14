import os
import sys
import glob
import itertools
import subprocess
import random
from itertools import repeat
from itertools import product
from multiprocessing import Pool

def dwgsim(arguments):
    error_rates = arguments[0]
    maf = arguments[1]
    coverage = arguments[2]
    genomes = arguments[3]
    seed = arguments[4]
    names = 'NeST_insilico'
    size = 20051 
    reads = int(size*coverage/500)
    dcmd = ['dwgsim', '-1', '250', '-2', '250', '-d', '1000', '-e', str(error_rates), '-F', str(maf), 
            '-E', str(error_rates), '-z', str(seed), '-c', '0', '-C', '-1', '-N', str(reads), '-I', '2', 
            genomes, '{0}{1}'.format(names, seed)]
    drun = subprocess.Popen(dcmd, shell=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    drun.wait()
    if drun.returncode != 0:
        print(' '.join(dcmd))
    else:
        return([names, genomes, coverage, reads, error_rates, maf, seed])

def cleanup():
    bfasts = glob.glob('*.bfast.fastq')
    vcfs = glob.glob('*.mutations.vcf')
    bwar1s = sorted(glob.glob('*.bwa.read1.fastq'))
    bwar2s = sorted(glob.glob('*.bwa.read2.fastq'))
    txts = glob.glob('*.mutations.txt')
    if not os.path.exists('./mutations'):
        os.mkdir('./mutations')
    for vcf, txt in zip(vcfs, txts):
        os.rename(vcf, './mutations/{0}'.format(vcf))
        os.rename(txt, './mutations/{0}'.format(txt))
    for bfast in bfasts:
        os.remove(bfast)
    for r1, r2 in zip(bwar1s, bwar2s):
        name = ''.join(r1.split('.')[:1])
        os.rename(r1, '{0}_r1.fastq'.format(name))
        os.rename(r2, '{0}_r2.fastq'.format(name))
        

def dgparallel():
    study = open('Plasmodium_insilico_study.tsv')
    error_rates = list()
    maf = list()
    coverage = list()
    reference = list()
    seed = list()
    for index, lines in enumerate(study):
        if index == 0 :
            continue
        lines = lines.strip().split('\t')
        error_rates.append(float(lines[5]))
        maf.append(float(lines[6]))
        coverage.append(float(lines[3]))
        reference.append(lines[2])
        seed.append(int(lines[-1]))
    pools = Pool(4)
    details = pools.map(dwgsim, zip(error_rates, maf, coverage, reference, seed))
    cleanup()

if __name__ == '__main__':
    dgparallel()

