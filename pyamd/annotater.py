import os
import sys
import csv
import vcf
from collections import namedtuple
from collections import OrderedDict
from operator import attrgetter
from pyamd.reader import Reader


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
            if rec.chrom == 'MT':
                chrom = '{0}{1}'.format(rec.chrom, rec.gene)
            else:
                chrom = rec.chrom
            try:
                coding_dict[chrom] += fasta_dict[rec.chrom][rec.start-1 : rec.stop]
            except KeyError:
                coding_dict[chrom] = fasta_dict[rec.chrom][rec.start-1 : rec.stop]
    
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

    def getBedOrder(self, bed_path):
        bed_order = list()
        reader = Reader()
        bed_reader = reader.readBed(bed_path)
        for bed_rec in bed_reader:
            if bed_rec.chrom not in bed_order:
                bed_order.append(bed_rec.chrom)
        return(bed_order)

    def getVcfOrder(self, vcf_path):
        vcf_order = list()
        reader = Reader()
        vcf_reader = reader.readVcf(vcf_path)
        for vcf_rec in vcf_reader:
            if vcf_rec.chrom not in vcf_order:
                vcf_order.append(vcf_rec.chrom)
        return(vcf_order)

    def getRevComp(self, fasta):
        rev_comp = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
        rev = ''
        for nuc in fasta:
            rev += rev_comp[nuc]
        return(rev)

    def getAA(self, codon):
        codontable = {
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
            'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
            'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
            }

        return(codontable[codon])

    def iterVcf(self, bed_path, vcf_path, fasta_path,name):
        reader = Reader()
        out_path = '{0}/variants_{1}.bed'.format(self.outdir,name)
        out_file = open(out_path, 'w')
        out_file.write('Chrom\tPos\tRef\tAlt\tExon\tRefCodon\tAltCodon\tCodonNumber\tRefAA\tAltAA\tCoverage\tQual\tAF\tPval\n')
        bed_reader = reader.readBed(bed_path)
        vcf_reader = reader.readVcf(vcf_path)
        contigs = reader.getVcfLength(vcf_path)
        coding_dict = self.getCodingFasta(fasta_path, bed_path)
        bed_ord = self.getBedOrder(bed_path)
        vcf_ord = self.getBedOrder(bed_path)    #Need to fix this
        if vcf_ord != bed_ord:
            print('Bed or VCF file not in order')
        bed_rec = next(bed_reader)
        vcf_rec = next(vcf_reader)
        mrna_len = 0
        bed_changed = 0
        while True:
            try:
                #Check if vcf record and bed record have the same chromosome
                if vcf_rec.chrom == bed_rec.chrom:
                    if vcf_rec.pos < bed_rec.start:
                        anno = "Intron"
                        try:
                            alfreq = self.getAlFreq(vcf_rec.info.DP4)
                        except AttributeError:
                            alfreq = self.getAlFreq(vcf_rec.format['Test'].AD)
                        pval = 100**(vcf_rec.qual/float(-10))
                        #Wrie intronic var to file  
                        out_file.write('{0}\t{1}\t{2}\t{3}\t{4}\tNA\tNA\tNA\tNA\tNA\t{5}\t{6}\t{7}\t{8}\n'.format(vcf_rec.chrom, 
                                        vcf_rec.pos, vcf_rec.ref, vcf_rec.alt, anno, vcf_rec.info.DP, vcf_rec.qual, alfreq, pval))
                        prev_vcf = vcf_rec.chrom
                        vcf_rec = next(vcf_reader)
                        bed_changed = 0
                
                    elif vcf_rec.pos in range(bed_rec.start, bed_rec.stop + 1):
                        anno = bed_rec.gene
                        try:
                            alfreq = self.getAlFreq(vcf_rec.info.DP4)
                        except AttributeError:
                            alfreq = self.getAlFreq(vcf_rec.format['Test'].AD)
                        pval = 10**(vcf_rec.qual/float(-10))
                        if vcf_rec.chrom == 'MT':
                            fasta = coding_dict['{0}{1}'.format(bed_rec.chrom, bed_rec.gene)]
                        else:
                            fasta = coding_dict[bed_rec.chrom]
                        if bed_rec.strand == '+':
                            codon_pos = vcf_rec.pos - bed_rec.start + mrna_len + 1
                            alt = vcf_rec.alt
                            if len(vcf_rec.ref) > 1 or len(vcf_rec.alt) > 1 :
                                codon = ['NA', 'NA', 'NA']
                                refAA = 'NA'
                                altAA = 'NA'
                            else:
                                codon = self.getCodon(codon_pos, fasta, alt)
                                refAA = self.getAA(codon[0])
                                altAA = self.getAA(codon[1])
                            out_file.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\n'.format(vcf_rec.chrom, 
                                            vcf_rec.pos, vcf_rec.ref, vcf_rec.alt, anno, codon[0], codon[1], codon[2], refAA, altAA, vcf_rec.info.DP,
                                            vcf_rec.qual, alfreq, pval))
                        elif bed_rec.strand == '-':
                            fasta = self.getRevComp(fasta[::-1])
                            rev_comp = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
                            alt = rev_comp[vcf_rec.alt]
                            codon_pos = bed_rec.stop - vcf_rec.pos + mrna_len + 1
                            if len(vcf_rec.ref) > 1 or len(vcf_rec.alt) > 1 :
                                codon = ['NA', 'NA', 'NA']
                                refAA = 'NA'
                                altAA = 'NA'
                            else:
                                codon = self.getCodon(codon_pos, fasta, alt)
                                refAA = self.getAA(codon[0])
                                altAA = self.getAA(codon[1])
                            out_file.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\n'.format(vcf_rec.chrom, 
                                            vcf_rec.pos, rev_comp[vcf_rec.ref], rev_comp[vcf_rec.alt], anno, codon[0], codon[1], codon[2], 
                                            refAA, altAA, vcf_rec.info.DP, vcf_rec.qual,alfreq, pval))
                        prev_vcf = vcf_rec.chrom
                        vcf_rec = next(vcf_reader)       
                        bed_changed = 0
                    
                    elif vcf_rec.pos > bed_rec.stop :
                        if bed_rec.chrom == 'MT':
                            mrna_len = 0
                        else:
                            mrna_len += bed_rec.stop - bed_rec.start +1
                        prev_bed = bed_rec.chrom
                        bed_rec = next(bed_reader)
                        bed_changed = 1 

                elif vcf_rec.chrom != bed_rec.chrom and bed_ord.index(bed_rec.chrom) < vcf_ord.index(vcf_rec.chrom):
                    mrna_len = 0
                    prev_bed = bed_rec.chrom
                    bed_rec = next(bed_reader)
                    bed_changed = 1
            
                elif vcf_rec.chrom != bed_rec.chrom and bed_ord.index(bed_rec.chrom) > vcf_ord.index(vcf_rec.chrom):
                    mrna_len = 0
                    anno = "Intron"
                    try:
                        alfreq = self.getAlFreq(vcf_rec.info.DP4)
                    except AttributeError:
                        alfreq = self.getAlFreq(vcf_rec.format['Test'].AD)
                    pval = 100**(vcf_rec.qual/float(-10))
                    #Wrie intronic var to file  
                    out_file.write('{0}\t{1}\t{2}\t{3}\t{4}\tNA\tNA\tNA\tNA\tNA\t{5}\t{6}\t{7}\t{8}\n'.format(vcf_rec.chrom, 
                                    vcf_rec.pos, vcf_rec.ref, vcf_rec.alt, anno, vcf_rec.info.DP, vcf_rec.qual, alfreq, pval))
                    prev_vcf = vcf_rec.chrom
                    vcf_rec = next(vcf_reader)
                    bed_changed = 0

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
