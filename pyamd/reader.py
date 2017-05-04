import re
import os
import sys
import csv
import vcf
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
            extract[lines.header] = lines.seq
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


    def getContigIDs(self, rec):
        '''Given a contig header record, return parsed contig object'''
        rec = re.split('##|=|<|,|>', rec)
        Contig = namedtuple('Contig', ['id', 'length'])
        contig = Contig(rec[4], int(rec[6]))
        return(contig)


    def getInfoIDs(self, rec):
        '''Given a INFO header records, return parsed INFO object'''
        rec = re.split('##|=|<|,|>|\"', rec)
        Info = namedtuple('Info', ['id', 'number', 'types', 'desc'])
        info = Info(rec[4], rec[6], rec[8], rec[11])
        return(info)

    def getFormatIDs(self, rec):
        '''Given a FORMAT header record, return parsed FORMAT object'''
        rec = re.split('##|=|<|,|>|\"', rec)
        Format = namedtuple('Format', ['id', 'number', 'types', 'desc'])
        formats = Format(rec[4], rec[6], rec[8], rec[11])
        return(formats)
    
    def getSamples(self, rec):
        '''Given the vcf header, return the number of samples and list of names'''
        rec = re.split('#|\t', rec)
        samples = rec[10:]
        number = len(samples)
        Sample = namedtuple('Sample', ['count', 'samples'])
        sample = Sample(number, samples)
        return(sample)

    def getRecInfo(self, info, info_dict):
        rec = re.split(';', info)
        info_val = {val.split('=')[0]: val.split('=')[1] for val in rec}
        info_keys = info_dict.keys()
        Info = namedtuple('Info', info_keys)
        values = list()
        for vals in info_dict:
            try:
                vals_id = vals
                vals_value = info_val[vals]
                if info_dict[vals_id].number == '4' and info_dict[vals_id].types == 'Integer':
                    values.append([int(allele) for allele in vals_value.split(',')])
                elif info_dict[vals_id].number == '4' and info_dict[vals_id].types == 'Float':
                    values.append([int(allele) for allele in vals_value.split(',')])
                elif info_dict[vals_id].number == 'A' and info_dict[vals_id].types == 'Integer':
                    values.append([int(val) for val in vals_value.split(',')])
                elif info_dict[vals_id].number == 'A' and info_dict[vals_id].types == 'Float':
                    values.append([float(val) for val in vals_value.split(',')])
                elif info_dict[vals_id].types == 'Integer':
                    values.append(int(vals_value))
                elif info_dict[vals_id].types == 'Float':
                    values.append(float(vals_value))
                elif info_dict[vals_id].types == 'Flag':
                    values.append(bool(vals_value))
                else:
                    values.append(str(vals_value))
            except KeyError:
                values.append(None)

        info = Info(*values)
        return(info)

    def getRecFormat(self, formats, format_dict, format_order):
        rec = re.split(':', formats)
        order = re.split(':', format_order)
        format_val = {keys : val for keys, val in zip(order,rec)}
        Format = namedtuple('Format', format_dict.keys())
        values = list()
        for vals in format_dict:
            try:
                vals_id = vals
                vals_value = format_val[vals]
                if format_dict[vals_id].number == 'G':
                    values.append(vals_value)
                elif format_dict[vals_id].number == '.':
                    values.append(vals_value)
                elif format_dict[vals_id].number == 'R':
                    values.append([int(allele) for allele in vals_value.split(',')])
                elif format_dict[vals_id].number == '4':
                    values.append([int(allele) for allele in vals_value.split(',')])
                elif format_dict[vals_id].types == 'Integer':
                    values.append(int(vals_value))
                elif format_dict[vals_id].types == 'Float':
                    values.append(float(vals_value))
                elif format_dict[vals_id].types == 'Flag':
                    values.append(bool(vals_value))
                else:
                    values.append(str(vals_value))
            except KeyError:
                values.append(None)

        formats = Format(*values)
        return(formats)

    def getVcfLength(self, vcf_path):
        contigs = OrderedDict()
        vcf_file = open(vcf_path, 'r')
        for lines in vcf_file:
            lines = lines.strip()
            if '##' == lines[:2]:
                if 'contig' == lines[2:8]:
                    rec_contig = self.getContigIDs(lines)
                    contigs[rec_contig.id] = rec_contig.length
            else:
                break
        return(contigs)

    def readVcf(self, vcf_path):
        vcf_file = open(vcf_path, 'r')
        Vcf = namedtuple('Vcf', ['chrom', 'pos', 'id', 'ref',
                'alt', 'qual', 'filter', 'info','format'])
        contigs = OrderedDict()
        info = OrderedDict()
        formats = OrderedDict()
        #samples = namedtuple()
        for lines in vcf_file:
            lines = lines.strip()
            print(lines)
            if '##' == lines[:2]:
                if 'contig' == lines[2:8]:
                    rec_contig = self.getContigIDs(lines)
                    contigs[rec_contig.id] = rec_contig.length
                elif 'INFO' == lines[2:6]:
                    rec_info = self.getInfoIDs(lines)
                    info[rec_info.id] = rec_info
                elif 'FORMAT' == lines[2:8]:
                    rec_format = self.getFormatIDs(lines)
                    formats[rec_format.id] = rec_format

            elif '#' == lines[0]:
                samples = self.getSamples(lines)
            
            elif '#' != lines[0]:
                var = lines.strip().split('\t')
                chrom = var[0]
                pos = int(var[1])
                ids = str(var[2])
                ref = var[3].split(',')[0]
                alt = var[4].split(',')[0]
                filters = var[6]
                if var[5] == '.':
                    qual = 0.0
                else:                    
                    qual = float(var[5])
                        
                infos = self.getRecInfo(var[7], info) 
                sample = OrderedDict()
                for sam, form in zip(samples.samples, var[9:]):
                    sample[sam] = self.getRecFormat(form, formats, var[8])
                record = Vcf(chrom, pos, ids, ref, alt, qual, filters, infos, sample)
                yield(record)

    def getAA(self, codon):
        codontable = {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
                      'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
                      'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
                      'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
                      'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
                      'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
                      'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
                      'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
                      'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
                      'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A'}
        change = records.alt + codon[1:]
        refAA = self.getAA(codon)
        altAA = self.getAA(change)
        record = Var(records.chrom, records.pos, records.ref, 
                records.alt, codon_pos, codon, change, refAA, altAA,
                records.depth, records.freq)
        yield record

            

if __name__ == '__main__':
    vcf_path = sys.argv[1]
    output = './'
    fastafile = Reader()
    reader = fastafile.readVcf(vcf_path)
    for lines in reader:
        print(lines)
