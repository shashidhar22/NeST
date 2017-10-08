import os
import sys
import csv
import time
import logging
import argparse
import subprocess

class Vcf:
    def __init__(self, vcf_path, bed_path, ref_path):
        self.bed_path = bed_path
        self.ref_path = ref_path
        return

    def getIns(self, sequence, loci_start, alternate, reference, offset):
     #   print(''.join(sequence))
     #   print(loci_start+1+offset, alternate)
        sequence.insert(loci_start+1+offset, alternate)
     #   print(sequence)
        loci_end = loci_start
        offset +=  1
        #print(loci_start, offset, alternate, reference, ''.join(sequence))
        return(sequence, loci_end, offset)

    def getDel(self, sequence, loci_start, alternate, reference, offset):
        corrected = list()
        deletion = iter(list(reference))

        for index, char in enumerate(sequence):
            if index >= (loci_start + offset) and index <= (loci_start + offset + len(reference) -1) :
                var = next(deletion)
                if var == char:
   #                 print(var)
                    continue
                else:
                    print('The nucleotides dont match {0}:{1} at {2}'.format(var, char, index))
            else:
                corrected.append(char)
        loci_end = loci_start + len(reference)
        offset += -len(reference)
        #print(loci_start, offset, alternate, reference, ''.join(sequence))
        return(corrected, loci_end, offset)

    def getSnp(self, sequence, loci_start, alternate, reference, offset):
  #      print(sequence, loci_start, alternate, reference, offset)
        if sequence[loci_start + offset] == reference:
  #          print('Reference base does match at {0}: {1} == {2}'.format(loci_start, reference, sequence[loci_start+offset]))
            sequence[loci_start + offset] = alternate
        else:
            print('Reference base doesnt match at {0}: {1} != {2}'.format(loci_start, reference, sequence[loci_start+offset]))
        loci_end = loci_start
        offset += 0
        #print(loci_start, offset, alternate, reference, ''.join(sequence))
        return(sequence, loci_end, offset)



   def rebuildRef(self, kestrel_path, reference, input_type, input_files, mindepth, ksize):
        out_fasta = '{0}/corrected.fa'.format(self.out_path)
        out_kestrel = '{0}/variants.tsv'.format(self.out_path)
        reader = Reader()
        kestrel = Kestrel(kestrel_path, reference, input_type, input_files,
                        self.out_path, mindepth, ksize)
        fasta_writer = FastaWriter(out_fasta)
        variants = kestrel.runKestrel()
        fasta = OrderedDict()
        for lines in reader.readFasta(reference):
            fasta[lines.header] = list(lines.seq)
        Fasta = namedtuple('Fasta',['header','seq'])
        record_index = 0
        loci_end = -1
        offset = 0
        gene = None
        for lines in variants:
            if lines.reference != gene:
                gene = lines.reference
                loci_end = -1
                offset = 0
            loci_start = lines.locus -1
            if lines.ref_depth != lines.var_depth:
 #               print('Not a homogenous var')
                continue
            elif lines.ref != '' and lines.alt == '' and loci_start > loci_end:
                fasta[gene], loci_end, offset = self.getDel(fasta[gene], loci_start, lines.alt, lines.ref, offset)
                #continue
            elif lines.ref == '' and lines.alt != '' and loci_start >= loci_end:
                fasta[gene], loci_end, offset = self.getIns(fasta[gene], loci_start, lines.alt, lines.ref, offset)
            elif lines.ref != '' and lines.alt != '' and loci_start > loci_end:
                fasta[gene], loci_end, offset = self.getSnp(fasta[gene], loci_start, lines.alt, lines.ref, offset)
#            print(gene, loci_start, offset, lines.ref, lines.alt, ''.join(fasta[gene]), lines.ref_depth, lines.var_depth)
        for lines in fasta:
            record = Fasta(lines, ''.join(fasta[lines]))
            fasta_writer.writeRec(record)

        return(out_fasta)

