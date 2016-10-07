import os
import sys
import csv
import vcf
from collections import namedtuple
from collections import OrderedDict
from operator import attrgetter
from reader import Reader


class Annotate:

    def __init__(self, outdir):
        self.outdir = os.path.abspath(outdir)
        return
    
    def getCodon(self, codon_pos, fasta, alt):
        '''Get the codon change for the variant.'''
        #The assumption made here is that the CDS has just on ORF
        codon = str()
        codon_number = int()
        if codon_pos % 3 ==0:
            codon = fasta[codon_pos-3: codon_pos]
            alt = codon[0:2] + alt
            codon_number = codon_pos / 3
        elif (codon_pos + 1) % 3 == 0:
            codon = fasta[codon_pos-2: codon_pos+1]
            alt = codon[0] + alt + codon[-1]
            codon_number = (codon_pos+1)/3
        elif (codon_pos +2) % 3 == 0:
            codon = fasta[codon_pos -1 : codon_pos +2]
            alt= alt + codon[1:]
            codon_number = (codon_pos+2)/3
        return(codon, alt, int(codon_number))

    def getCodingFasta(self, fasta_path, bed_path):
        '''Given a bed file with exons and reference fasta.
        Return the CDS for each gene.'''
        reader = Reader()
        bed_reader = reader.readBed(bed_path)
        fasta_dict = reader.extractFasta(fasta_path)
        coding_dict = OrderedDict()
        for rec in bed_reader:
            try:
                coding_dict[rec.chrom] += fasta_dict[rec.chrom][rec.start-1 : rec.stop]
            except KeyError:
                coding_dict[rec.chrom] = fasta_dict[rec.chrom][rec.start-1 : rec.stop]
    
        return(coding_dict)

    def getAlFreq(self, depth):
        '''Calculate allele frequency based on allelic depth'''
        total = sum(depth)
        if len(depth) == 4:
            ref = sum(depth[:2])
            alt = sum(depth[2:])
            alfreq = alt/float(total)
        else:
            ref = depth[0]
            alt = depth[1]
            alfreq = alt/float(total)
        return(alfreq) 

    def iterVcf(self, bed_path, vcf_path, fasta_path,name):
        reader = Reader()
        out_path = '{0}/variants_{1}.bed'.format(self.outdir,name)
        out_file = open(out_path, 'w')
        out_file.write('Chrom\tPos\tRef\tAlt\tExon\tRefCodon\tAltCodon\tCodonNumber\tAF\tPval\n')
        bed_reader = reader.readBed(bed_path)
        vcf_reader = reader.readVcf(vcf_path)
        coding_dict = self.getCodingFasta(fasta_path, bed_path)
        bed_rec = next(bed_reader)
        vcf_rec = next(vcf_reader)
        mrna_len = 0
        while True:
            try:
                #Change bed record if var chromosome does not match
                if vcf_rec.chrom != bed_rec.chrom:
                    mrna_len = 0
                    bed_rec = next(bed_reader)
            
                elif vcf_rec.chrom != bed_rec.chrom: 
                    mrna_len = 0 
                    vcf_rec = next(vcf_reader)

                #Change bed record if var position is beyond bed stop
                elif vcf_rec.pos > bed_rec.stop:
                    mrna_len += bed_rec.stop - bed_rec.start +1
                    bed_rec = next(bed_reader)
                

                #If var position is before bed start, annotate var as Intronic variants            
                elif vcf_rec.pos < bed_rec.start:
                    anno = "Intron"
                    try:
                        alfreq = self.getAlFreq(vcf_rec.info.DP4)
                    except AttributeError:
                        alfreq = self.getAlFreq(vcf_rec.format['Test'].AD)
                    pval = 10**(vcf_rec.qual/float(-10))
                    #Write intronic var to file
                    out_file.write('{0}\t{1}\t{2}\t{3}\t{4}\tNA\tNA\tNA\t{5}\t{6}\n'.format(vcf_rec.chrom, 
                                    vcf_rec.pos, vcf_rec.ref, vcf_rec.alt, anno, alfreq, pval))
                    vcf_rec = next(vcf_reader)
                #If var in bed range, annotate with exon number and get codon change
                elif vcf_rec.pos in range(bed_rec.start, bed_rec.stop + 1) :
                    anno = bed_rec.gene
                    fasta = coding_dict[vcf_rec.chrom]
                    codon_pos = vcf_rec.pos - bed_rec.start + mrna_len +1
                    codon = self.getCodon(codon_pos, fasta, vcf_rec.alt)
                    try:
                        alfreq = self.getAlFreq(vcf_rec.info.DP4)
                    except AttributeError:
                        alfreq = self.getAlFreq(vcf_rec.format['Test'].AD)
                    pval = 10**(vcf_rec.qual/float(-10))
                    out_file.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n'.format(vcf_rec.chrom, 
                                    vcf_rec.pos, vcf_rec.ref, vcf_rec.alt, anno, codon[0], codon[1], codon[2], 
                                    alfreq, pval))
                    vcf_rec = next(vcf_reader)
            except StopIteration:
                break
        out_file.close()
        return


if __name__ == '__main__':
    bed_path = sys.argv[1]
    fasta_path = sys.argv[2]
    vcf_path = sys.argv[3]
    out_path = os.path.abspath(sys.argv[4])
    annotate = Annotate(out_path)
    name = sys.argv[5]
    annotate.iterVcf(bed_path, vcf_path, fasta_path, name)
