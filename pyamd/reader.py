import os
import sys
import csv
from collections import namedtuple
from collections import OrderedDict
from operator import attrgetter

class Reader:

    def __init__(self):
        return

    def readFasta(self, files):
        '''Given a fasta file, read the iterate through the file and yield 
        each record, along with the header information, custom fasta id,
        length of sequence. If minflag is set to True, also calculate the 
        minimizer for the sequence.'''
        fasta_handle = open(files)
        #Read whole file
        #FIXME : Reading the whole file to memory can go out of hand.
        #TODO: Implement a better iterator ; logged 07/22/2016
        fasta_content = fasta_handle.read().split('>')[1:]
        fastaid = 0
        #FIXME : Using two named tuples to account for minimizer may not be the best approach.
        #TODO: Acount for minimizer flag in a different manner ; logged 07/22/2016
        Fasta = namedtuple('Fasta',['header','seq','fid','length'])
        for lines in fasta_content:
            header = lines.split('\n')[0].split(' ')[0]
            sequence = ''.join(lines.split('\n')[1:]).upper()
            length = len(sequence)
            fastaid = lines.split('\n')[0].split('_')[0]
            record = Fasta(header, sequence, fastaid, length)
            yield record

    def filterFasta(self, contig_list, fasta):
        fastafile = self.readFasta(fasta)
        for lines in fastafile:
            if lines.header in contig_list:
                continue
            yield lines

    def readBed(self, bed_path):
        '''Given a bed file extract chromsome, start, stop, genename, score 
        and strand information. Returns a namedtuple for each record in bed 
        file.'''
        bed_handle = open(bed_path)
        Bed = namedtuple('Bed', ['chrom','start','stop','gene','score', 'strand'])
        for lines in bed_handle:
            line = lines.strip()
            print(line.split('\t'))
            chrom = line.split('\t')[0]
            start = int(line.split('\t')[1])
            stop = int(line.split('\t')[2])
            gene = line.split('\t')[3].split(',')[0]
            score = lines.split('\t')[4]
            strand = line.split('\t')[5]
            record = Bed(chrom, start, stop, gene, score, strand)
            yield record

    def extractFasta(self, fasta_file):
        '''Given a reference fasta extract seqeunces. Returns an ordered 
        dictionary, with gene names as key and sequence as value.'''
        #FIXME: Storing sequences in dictionaries can get out of hand.
        #TODO: Implement lazy loading of extracted sequences ; logged 07/22/2016.
        fasta = self.readFasta(fasta_file)
        Fasta = namedtuple('Fasta',['header','seq','fid','length'])
        extract = OrderedDict()
        for lines in fasta:
            extract[lines.header] = lines
        return(extract)

    def getGenes(self, fasta_file):
        '''Given a bed file, get all the genes listed in the file. Returns
        a set containing all the gene names.'''
        fasta = self.readFasta(fasta_file)
        genes = set()
        for lines in fasta:
            genes.add(lines.header)
        return(genes)
            
    def writeExtract(self, fasta_file):
        '''Given extracted sequences, write the sequences to a file'''
        #TODO: Move to a write class; logged 07/22/2016
        ext = self.extractFasta()
        output = open('{0}/extracted.fa'.format(self.outpath), 'w')
        for lines in ext:
            output.write('>{0}\n{1}\n'.format(ext[lines].header, ext[lines].seq))
        output.close()
        return

    def readVcf(self, vcf_path):
        vcf_file = open(vcf_path, 'rt')
        Vcf = namedtuple('Vcf', ['chrom', 'pos', 'id', 'ref',
                'alt', 'qual', 'filter', 'info', 'sample'])
        for lines in vcf_file:
            if lines[0] == '#':
                continue
            var = lines.strip().split('\t')
            chrom = var[0]
            pos = int(var[1])
            ids = var[2]
            ref = var[3]
            alt = var[4]
            qual = float(var[5])
            filters = var[6]
            info = var[7]
            sample = var[8]
            record = Vcf(chrom, pos, ids, ref, alt, qual, filters,
                        info, sample)
            yield record

    def extractVars(self, vcf_path, fasta_path):
        vcf = self.readVcf(vcf_path)
        Var = namedtuple('Var', ['chrom', 'pos', 'ref', 'alt', 'codonPos','refCodon', 'altCodon'])
        fasta = {rec.header : rec.seq for rec in self.readFasta(fasta_path)}
        for records in vcf:
            if (records.pos) % 3 == 0:
                codon_pos = (records.pos)/3
                codon = fasta[records.chrom][(records.pos-3):records.pos]#[(records.pos-2):records.pos+1]
                print(records, codon)
                change = codon[:2] + records.alt 
                record = Var(records.chrom, records.pos, records.ref, records.alt, codon_pos, codon, change)
                yield record
            elif (records.pos +1 ) % 3 == 0:
                codon_pos = (records.pos +1)/3
                codon = fasta[records.chrom][(records.pos-2):(records.pos+1)]#[records.pos-1:records.pos+2]
                print(codon)
                change = codon[0] + records.alt + codon[2]
                record = Var(records.chrom, records.pos, records.ref, records.alt, codon_pos, codon, change)
                yield record
            elif (records.pos +2) % 3 == 0:
                codon_pos = (records.pos +2)/3
                codon = fasta[records.chrom][(records.pos-1):(records.pos+2)]#[records.pos:records.pos+3]
                print(codon)
                change = records.alt + codon[1:]
                record = Var(records.chrom, records.pos, records.ref, records.alt, codon_pos,codon, change)
                yield record

            

if __name__ == '__main__':
    fasta = sys.argv[1]
    bed = ''
    output = './'
    fastafile = Reader(bed, fasta, output)
    reader = fastafile.readFasta()
    for lines in sorted(reader, key=attrgetter('minimizer')):
        print(lines.minimizer)


