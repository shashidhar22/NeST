import os
import sys
import vcf
import glob
import math
import pandas as pd
import pysam
from parsers import Fasta
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

class Summary:

    def __init__(self, fasta, bed, voi, out_path):
        self.fasta = fasta
        self.bed = bed
        self.voi = voi
        self.out_path = out_path

    def createTable(self):
        voi = pd.read_excel('./Reportable_SNPs_Report_v2.xlsx', sheetname=1, parse_cols="B:C")
        voi['Ref'] = voi['SNP'].str.split('\d+', expand=True)[0]
        voi['Alt'] = voi['SNP'].str.split('\d+', expand=True)[1]
        voi['Pos'] = voi['SNP'].str.split('[a-zA-Z]',expand=True)[1]
        voi = voi[['Gene','Pos','Ref','Alt']]
        voi[['Pos']] = voi[['Pos']].apply(pd.to_numeric)
        voi.set_index(['Gene','Pos'], inplace=True)
        # Read fasta file
        reader = Fasta('../../ref/mdr.fa')
        #Read bed file
        bed = pd.read_table('../../ref/mdr.bed', header=None, sep='\t', names=['Gene','Start','Stop','Exon','Info','Strand'])
        #Spliting in codons
        bed_iter = iter(bed.iterrows())
        fasta = reader.read()
        codon_dict = dict()
        codon = list()
        codon_table = list()
        for records in fasta:
            seq = iter(records.seq)
            cdna = []
            for pos, base in enumerate(seq, 1):
                for exons in bed.iterrows():
                    if exons[1]['Gene'] == records.header:
                        if pos in range(exons[1]['Start'], exons[1]['Stop']+1):
                            cdna.append((base, pos, exons[1]['Strand']))
                        else:
                            continue
                    else:
                        continue
            aa = []
            codon = 0
            sequence = cdna #iter(cdna)
            for bases in range(0, len(cdna),3):
                codon += 1
                codon_list = cdna[bases:bases+3] if cdna[bases][2] == '+' else cdna[bases:bases+3][::-1]
                codon_seq = ''.join([val[0] for val in codon_list])
                for codon_pos, nuc in enumerate(codon_list):
                    codon_table.append((records.header, (records.header, nuc[1], nuc[0], codon, codon_seq, getAA(codon_seq), codon_pos)))

        codon_data = pd.DataFrame.from_items(codon_table, orient='index', columns=['Gene', 'NucPos', 'Nuc', 'Pos', 'Codon', 'AA', 'CodonPos'])
        codon_data[['Pos']] = codon_data[['Pos']].apply(pd.to_numeric)
        codon_data.set_index(['Gene', 'Pos'], inplace=True)
        codons_of_int = voi.merge(codon_data, left_index=True, right_index=True, how='inner' )
        codons_of_int['AltCodon'] = codons_of_int['Alt'].apply(getCodon)
        return(codons_of_int)

    def getBamStat(self, bamfile, chrom, pos, codon_pos):
        bamfile = pysam.AlignmentFile(bamfile, 'rb')
        codon_coverage = 0
        if codon_pos == 0:
            avg_codon_coverage = bamfile.count(chrom, pos-1, pos+2)
        elif codon_pos == 1:
            avg_codon_coverage = bamfile.count_coverage(chrom, pos-2, pos+1)
        else:
            avg_codon_coverage = samfile.count_coverage(chrom, pos-3, pos)
        return(avg_codon_coverage)

    def getVarOfInt(self, vcf_path):
        sample_dirs = glob.glob('{0}/*'.format(out_path))
        for samples in sample_dirs:
            sample_folder = os.path.basename(samples)
            sample_bam = '{0}/output_sorted_RG.bam'
            sample_vcf = vcf.Reader(filename='{0}/{0}_merged_annotated.vcf'.format(samples))
            sample_name = sample_vcf.samples[0]
            sample_result = []
            sample_vars = {'{0}_{1}_{2}_{3}'.format(var.CHROM, var.POS, var.INFO['RefAA'][0], var.INFO['AltAA'][0]): var.INFO['DP'] for var in sample_vcf}
            for val in codons_of_int.iterrows():
                if '{0}_{1}_{2}_{3}'.format(val[0][0], val[1]['NucPos'], val[1]['Ref'],val[1]['Alt']) in sample_vars.keys():
                    depth = sample_vars['{0}_{1}_{2}_{3}'.format(val[0][0], val[1]['NucPos'], val[1]['Ref'],val[1]['Alt'])]
                    sample_result.append(depth)
                else:
                    depth = -1 * self.getBamStat(sample_bam, val[0][0], val[1]['NucPos'], val[1]['CodonPos'])
                    sample_result.append(depth)
            sample_series = pd.Series(sample_result, index=codons_of_int.index)
            codons_of_int[sample_name] = sample_series
        sns.set()
        sns.set_context("notebook")
        plt.figure(figsize=(24,30))
        heatmap = sns.heatmap(codons_of_int.iloc[:,9:],linewidths=.5, square=False)
        fig = heatmap.get_figure()
        fig.savefig('{0}/heatmap.png'.format(self.out_path))


if __name__ == '__main__':
	fasta_path = sys.argv[1]
    bed_path = sys.argv[2]
    voi_path = sys.argv[3]
	out_path = sys.argv[4]
	summarizer = Summary(fasta_path, bed_path, voi_path, out_path)
	summarizer.getVarOfInt(out_path)
    
