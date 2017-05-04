
# coding: utf-8

# In[1]:

import os
import sys
import vcf
import glob
import pandas as pd


# ## Custom parsers

# ### Fasta parser

# In[2]:

import os
import re
import sys
import time
import glob
import numpy
import argparse
from operator import attrgetter
from collections import namedtuple
from collections import OrderedDict

class Fasta:

	def __init__(self, fasta_path):
		self.fasta_path = fasta_path
		return

	def peek(self, fasta_handle):
		curr_pos = fasta_handle.tell()
		curr_line = fasta_handle.readline()
		fasta_handle.seek(curr_pos)
		return(curr_line)

	def read(self):
		'''Given a fasta file, read will iterate through the file and yield each record, alonf with the header information, custom fasta id, length of sequence.'''
		fasta_handle = open(self.fasta_path)

		#Read the file one line at a time and process a chunk of text until the next header is found. At this point, process the chunk of text as one sequence and create iterator.
		next_line = ''
		line_number = 0
		header = ''
		sequence = ''
		header_found = False
		fid = 0
		while True:
			fasta = namedtuple('fastaRec', ['header', 'seq', 'fid', 'length'])
			next_line = self.peek(fasta_handle)
			try:
					if next_line[0] == '>':
						header = fasta_handle.readline().strip()[1:]
					sequence = ''
					fid += 1
					header_found = True
					line_number += 1
					while True:
						try:
								if self.peek(fasta_handle)[0] != '>' and self.peek(fasta_handle) != ' \n': 
									sequence += fasta_handle.readline().strip()
									header_found = False
									line_number += 1

								elif (self.peek(fasta_handle) == ' \n' or self.peek(fasta_handle)[0] == '>') and header_found:
									line_number += 1
									raise SyntaxError('Sequence missing for header : {0} at line {1}'.format(header, line_number))
									sys.exit()
								elif self.peek(fasta_handle)[0] == '>' and not header_found:
									break					
						except IndexError:
								break								
					length = len(sequence)
					record = fasta(header, sequence, fid, length)
					yield record
			except IndexError:
				break					

	def write(self, out_path, reader_obj, wrapping = 0):
		'''Given a fasta object and a output path, will write out a fasta file.'''

		fasta_handle = open(out_path, 'w')
		seq = ''
		for sequences in reader_obj:
			if wrapping == 0:
				wrapping = len(sequences.seq)
			seq = re.findall('.{{1,{0}}}'.format(wrapping), sequences.seq)
			fasta_handle.write('>{0}|{1}|{2}\n'.format(sequences.header, sequences.length, len(''.join(seq))))
			for record in seq:
				fasta_handle.write('{0}\n'.format(record))
		fasta_handle.close()
		return					


# ### Codon table

# In[3]:

def getAA(codon):
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
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'
        }
    return(codontable[codon])


# ## Format variants of interest data

# In[4]:

voi = pd.read_excel('./Reportable_SNPs_Report_v2.xlsx', sheetname=1, parse_cols="B:C")
voi.head()


# In[5]:

voi['Ref'] = voi['SNP'].str.split('\d+', expand=True)[0]
voi['Alt'] = voi['SNP'].str.split('\d+', expand=True)[1]
voi['Pos'] = voi['SNP'].str.split('[a-zA-Z]',expand=True)[1]
voi = voi[['Gene','Pos','Ref','Alt']]
voi.set_index(['Gene','Pos'], inplace=True)
voi.head()


# ## Create codon table for genes

# In[11]:

# Read fasta file
reader = Fasta('../../ref/mdr.fa')
fasta = reader.read()

#Read bed file
bed = pd.read_table('../../ref/mdr.bed', header=None, sep='\t', names=['Gene','Start','Stop','Exon','Info','Strand'])
bed.head()


# In[12]:

#Spliting in codons
bed_iter = iter(bed.iterrows())
bed_rec = next(bed_iter)
fasta_rec = next(fasta)
fasta_pos = 0
fasta_seq = iter(fasta_rec)
codon_pos = 0
stop = bed_rec[1]['Stop'] - (bed_rec[1]['Stop'] %3)
carry = bed_rec[1]['Stop'] % 3
stop


# In[13]:

bed_rec


# In[14]:

fasta_rec


# In[16]:

while True:
    if fasta_rec.header == bed_rec[1]['Gene']:
        try:
            if fasta_pos < bed_rec[1]['Start']:
                print('In If')
                next(fasta_seq)
                fasta_pos += 1 # bed_rec[1]['Start']
            elif fasta_pos >= bed_rec[1]['Start'] and fasta_pos < stop:
                print('In elif 1')
                tmp = fasta_pos
                codon =  next(fasta_seq)+next(fasta_seq)+next(fasta_seq) #fasta_seq[(fasta_pos-1): (fasta_pos+2)] #next(fasta_seq)+next(fasta_seq)+next(fasta_seq)
                fasta_pos += 3
                codon_pos += 1
                print(bed_rec[1]['Gene'], codon, tmp, fasta_pos, codon_pos, getAA(codon))
            elif fasta_pos >= stop:
                print('In elif 2')
                codon = fasta_seq[stop: bed_rec[1]['Stop']]
                print(codon)
                bed_rec = next(bed_iter)
                print(bed_rec)
                stop = bed_rec[1]['Stop'] - (bed_rec[1]['Stop'] %3)
                print(stop)
                carry = 3 - (bed_rec[1]['Stop'] % 3)
                print(carry)
                fasta_pos = bed_rec[1]['Start']
                print(fasta_pos)
                while carry > 0:
                    codon += next(fasta_seq)
                    carry -= 1
                codon += fasta_seq[(fasta_pos-1): (fasta_pos + carry -1) ]
                print(codon)
                fasta_pos += carry
                print(fasta_pos)
                codon_pos += 1
                print(codon_pos)
                print(bed_rec[1]['Gene'], codon, tmp, fasta_pos, codon_pos, getAA(codon))
        except KeyError:
            fasta_rec = next(fasta)
            fasta_pos = 0 
            fasta_seq = iter(fasta_rec.seq)
            codon_pos = 0
    else


# In[ ]:



