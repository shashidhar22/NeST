import re
import os
import sys
import csv
import gzip
import time
import math
import logging
import numpy as np
from pprint import pprint
from collections import namedtuple
from collections import OrderedDict
from pyamd.parsers.fasta import Fasta


class Bed:

    def __init__(self, bed_path):
        self.bed_path = os.path.abspath(bed_path)
        self.getUid()

    def read(self):
        Bed = namedtuple('Bed', ['chrom', 'start', 'stop', 'name', 'score',
                                 'strand', 'thickStart', 'thickEnd', 'itemRGB',
                                 'blockCount', 'blockSizes', 'blockStarts',
                                 'length'])
        bed_handle = open(self.bed_path)
        for lines in bed_handle:
            line = lines.strip().split('\t')
            chrom = line[0]
            start = int(line[1])
            stop = int(line[2])
            length = stop - start + 1
            name = line[3]
            if line[4] == '.':
                score = 0
            else:
                score = int(line[4])
            strand = line[5]
            thickStart = int(line[6])
            thickEnd = int(line[7])
            itemRGB = None
            blockCount = int(line[9])
            blockSizes = [int(val) for val in line[10].split(',') if val != '']
            blockStarts = [int(val) for val in line[11].split(',') if val != '']
            if (len(blockSizes) != blockCount or
                len(blockStarts) != blockCount):
                print('Error in BED format, please check the block details')
                print(line)
                sys.exit()
            record = Bed(chrom, start, stop, name, score, strand, thickStart,
                        thickEnd, itemRGB, blockCount, blockSizes, blockStarts,
                        length)
            yield(record)

    def getUid(self):
        self.uids = OrderedDict()
        index = 1
        for records in self.read():
            if records.chrom not in self.uids:
                length = records.length
                order = 10**(int(math.log10(length)) + 1) * index
                index += 1
                self.uids[records.chrom] = order
        return

    def getExonTable(self):
        reader = self.read()
        Annotations = namedtuple('Table', ['chrom', 'gene', 'exon', 'start',
                                        'stop', 'strand', 'uidStart',
                                        'uidStop', 'codonStart', 'codonStop',
                                        'length', 'blockSizes'])
        for records in reader:
            exon_count = 1
            for estart, length in zip(records.blockStarts, records.blockSizes):
                chrom = records.chrom
                gene = records.name
                exon = 'exon{0}'.format(exon_count)
                start = estart
                stop = estart + length
                strand = records.strand
                uidStart = self.uids[chrom] + start
                uidStop = self.uids[chrom] + stop
                codonStart = records.thickStart
                codonStop  = records.thickEnd
                blockSize = records.blockSizes
                annotation = Annotations(chrom, gene, exon, start, stop, strand,
                                uidStart, uidStop, codonStart, codonStop,
                                length, blockSize)
                exon_count += 1
                yield annotation

    def checkOrder(self, fasta_path):
        fasta_file = Fasta(fasta_path)
        fasta_reader = fasta_file.read()
        fasta_dict = OrderedDict()
        for record in fasta_reader:
            fasta_dict[record.header] = record.fid

        for fasta_header, bed_header in zip(fasta_dict, self.uids):
            if fasta_dict[fasta_header] == self.uids[bed_header]:
                return(True)
            else:
                return(False)

    def getCodingFasta(self, fasta_path):
        fasta_file = Fasta(fasta_path)
        fasta_reader = fasta_file.read()
        exon_reader = self.getExonTable()
        if not self.checkOrder(fasta_path):
            print('Bed and fasta contigs are not in the same order')
            sys.exit()
        coding_seq = namedtuple('Fasta', ['chrom', 'seq', 'fid', 'length'])
        coding_dict = OrderedDict()
        fasta_rec = next(fasta_reader)
        exon_rec = next(exon_reader)
        while True:
            try:
                if exon_rec.uidStart in range(fasta_rec.fid,
                    fasta_rec.fid + fasta_rec.length + 1):
                    if exon_rec.gene not in coding_dict:
                        chrom = exon_rec.gene
                        seq = fasta_rec.seq[exon_rec.start-1: exon_rec.stop]
                        fid = exon_rec.uidStart
                        length = len(seq)
                        strand = exon_rec.strand
                        coding_dict[chrom] = [chrom, seq, fid, length, strand]
                    else:
                        chrom = exon_rec.gene
                        coding_dict[chrom][1] += fasta_rec.seq[exon_rec.start-1:
                                                            exon_rec.stop]
                        coding_dict[chrom][3] += len(coding_dict[chrom][1])
                    exon_rec = next(exon_reader)
                else:
                    fasta_rec = next(fasta_reader)
            except StopIteration:
                break

        for genes in coding_dict:
            chrom = coding_dict[genes][0]
            if coding_dict[genes][4] == '-':
                seq = self.getRevComp(coding_dict[genes][1])
            elif coding_dict[genes][4] == '+':
                seq = coding_dict[genes][1]
            fid = coding_dict[genes][2]
            length = coding_dict[genes][3]
            record = coding_seq(chrom, seq, fid, length)
            yield(record)

    def getRevComp(self, fasta):
#        logger.debug('Reverse complimenting fasta sequence')
        rev_comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        rev = ''
        for nuc in fasta:
            rev += rev_comp[nuc]
        return(rev)
