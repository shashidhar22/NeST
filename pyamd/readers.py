import os
import sys
import csv
import gzip
import time
import logging
import numpy as np
from collections import namedtuple


class Fastq:

    def __init__(self, fastqfile, outdir, phred):
        self.fastq = fastqfile
        self.outdir = outdir
        self.phred = phred
        FORMAT = '%(asctime)-15s : %(levelname)-8s :  %(message)s'
        logging.basicConfig(format=FORMAT)
        self.phreddict = self.phredmap()


    def phredmap(self):
        phreddict = dict()
        if self.phred == "phred33":
            for asciis, quals in zip(range(33,126), range(0,92)):
                phreddict[asciis] = quals
        return(phreddict)

    def formatChecker(self, header, seq, sheader, quals, line, phred):
        if header[0] != "@":
            return(1)
        if sheader == "" or sheader[0] != "+":
            return(2)
        if len(seq) != len(quals):
            return(3)
        if seq == '':
            return(4)

    def qualmasker(self, seq, quals):
        maskedseq = list()
        for base, qual in zip(seq, quals):
            if qual < 20:
                maskedseq.append('N')
            else:
                maskedseq.append(base)
        return(maskedseq, quals)

    def read(self):
        if '.gz' in self.fastq or '.fastqgz' in self.fastq :
            fqhandle = gzip.open(self.fastq, 'rb')
        else:
            fqhandle = open(self.fastq, 'r')
        Record = namedtuple('Record',['header', 'sheader', 'seq', 'quals'])
        lineno = 0
        while True:
            try:
                #Get fastq record
                if '.gz' in self.fastq or '.fastqgz' in self.fastq:
                    header = str(next(fqhandle),'utf-8').strip()
                    seq = [base for base in str(next(fqhandle),'utf-8').strip()]
                    sheader = str(next(fqhandle), 'utf-8').strip()
                    quals = str(next(fqhandle), 'utf-8').strip()
                else:
                    header = next(fqhandle).strip()
                    seq = [base for base in next(fqhandle).strip()]
                    sheader = next(fqhandle).strip()
                    quals = next(fqhandle).strip()
                #quals = [self.phreddict[int(ord(str(qual)))] for qual in next(fqhandle).strip()]
                lineno += 1

                #Check if record adheres to fastq format
                check = self.formatChecker(header, seq, sheader, quals, lineno, self.phred)
                if check == 1:
                    logging.error('Invalid header in fastq read ; record number : {0}'.format(lineno))
                    raise NameError('Invalid header in fastq read; record number : {0}'.format(lineno))
                if check == 2:
                    logging.error('Invalid secondary header in fastq read ; record number : {1}'.format(lineno))
                    raise NameError('Invalid secondary header in fastq read; record number : {0}'.format(lineno))
                if check == 3:
                    logging.error('Sequence and quality strings of unequl length in fastq read; record number : {0}'.format(lineno))
                    raise NameError('Sequence and quality strings of unequal length in fastq read; record number : {0}'.format(lineno))
                if check == 4:
                    logging.error('Sequence data missing; record number : {0}'.format(lineno))
                    raise NameError('Sequence data missing; record number : {0}'.format(lineno))

                #Optionally mask low quality reads
                #seq, quals = self.qualmasker(seq, quals)

                #Return record
                record = Record(header, sheader, seq, quals)
                yield(record)
            except StopIteration:
                logging.info('End of file')
                break
            except NameError:
                break

        return

class Bed:

    def __init__(self, bed_path):
        self.bed_path = os.path.abspath(bed_path)

    def read(self):
        Bed = namedtuple('Bed', ['chrom', 'start', 'stop', 'gene', 'score', 'strand'])
        bed_handle = open(self.bed_path)
        for lines in bed_handle:
            line = lines.strip().split('\t')
            chrom = line[0]
            start = int(line[1])
            stop = int(line[2])
            gene = line[3]
            if line[4] == '.':
                score = 0
            else:
                score = int(line[4])
            strand = line[5]
            record = Bed(chrom, start, stop, gene, score, strand)
            yield(record)
