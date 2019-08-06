import os
import re
import sys
import time
import glob
import math
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
		fid = 1
		while True:
			fasta = namedtuple('fastaRec', ['header', 'seq', 'fid', 'length'])
			next_line = self.peek(fasta_handle)
			try:
					if next_line[0] == '>':
						header = fasta_handle.readline().strip()[1:].split('|')[0].strip()
					sequence = ''
					header_found = True
					line_number += 1
					while True:
						try:
								if self.peek(fasta_handle)[0] != '>' and self.peek(fasta_handle) != ' \n':
									sequence += fasta_handle.readline().strip().upper()
									header_found = False
									line_number += 1

								elif (self.peek(fasta_handle) == ' \n' or self.peek(fasta_handle)[0] == '>') and header_found:
									line_number += 1
									raise SyntaxError('Sequence missing for header : {0} at line {1}'.format(header, line_number))
								elif self.peek(fasta_handle)[0] == '>' and not header_found:
									break
						except IndexError:
								break
					length = len(sequence)
					order = 10**(int(math.log10(length)) + 1) * fid
					fid += 1
					record = fasta(header, sequence, order, length)
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

if __name__ == '__main__':
	fasta_path = sys.argv[1]
	out_path = sys.argv[2]
	wrapping = int(sys.argv[3])
	fasta_reader = Fasta(fasta_path)
	fasta_reader.write(out_path, fasta_reader.read(), wrapping)
