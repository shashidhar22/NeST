import os
import sys
import vcf
import glob
import math
import pandas as pd
import pysam
from pyamd.parsers import Fasta
from pyamd.annotater import Annotate
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt

import seaborn as sns

class Summary:

    def __init__(self, fasta, bed, voi, out_path):
        self.fasta = fasta
        self.bed = bed
        self.voi = voi
        self.out_path = out_path

    def getVarStats(self, vcf_file):
        vcf_file = vcf.Reader(filename=vcf_file)
        total = 0
        exonic = 0
        intronic = 0
        verfied = 0
        syn = 0
        nsyn = 0
        trans = 0
        tranv = 0
        trasition = ['AG', 'GA', 'CT', 'TC']
        transversion = ['AC', 'AT', 'CA', 'CG', 'GC', 'GT', 'TA', 'TG']
        for variant in vcf_file:
            total += 1
            if variant.INFO['ExonNumber'][0] == 'Intron':
                intronic += 1
            else:
                exonic += 1
                if variant.INFO['Found'][0] == '2':
                    verfied += 1
                if variant.INFO['RefAA'] == variant.INFO['AltAA']:
                    syn += 1
                else:
                    nsyn += 1
                if '{0}{1}'.format(variant.REF, str(variant.ALT[0])) in trasition:
                    trans += 1
                else:
                    tranv += 1
        return(total, verfied, exonic, intronic, syn, nsyn, trans, tranv)

    def createTable(self):
        voi = pd.read_excel(self.voi, sheetname=1, parse_cols="B:C")
        voi['Ref'] = voi['SNP'].str.split('\d+', expand=True)[0]
        voi['Alt'] = voi['SNP'].str.split('\d+', expand=True)[1]
        voi['Pos'] = voi['SNP'].str.split('[a-zA-Z]',expand=True)[1]
        voi = voi[['Gene','Pos','Ref','Alt']]
        voi[['Pos']] = voi[['Pos']].apply(pd.to_numeric)
        voi.set_index(['Gene','Pos'], inplace=True)
        # Read fasta file
        annotater = Annotate(self.out_path)
        reader = Fasta(self.fasta)
        #Read bed file
        bed = pd.read_table(self.bed, header=None, sep='\t', names=['Gene','Start','Stop','Exon','Info','Strand'])
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
                    codon_table.append((records.header, (records.header, nuc[1], nuc[0], codon, codon_seq,
                                                         annotater.getAA(codon_seq), codon_pos)))

        codon_data = pd.DataFrame.from_items(codon_table, orient='index', columns=['Gene', 'NucPos', 'Nuc', 'Pos', 'Codon', 'AA', 'CodonPos'])
        codon_data[['Pos']] = codon_data[['Pos']].apply(pd.to_numeric)
        codon_data.set_index(['Gene', 'Pos'], inplace=True)
        codons_of_int = voi.merge(codon_data, left_index=True, right_index=True, how='inner' )
        codons_of_int['AltCodon'] = codons_of_int['Alt'].apply(annotater.aaToCodon)
        return(codons_of_int)

    def getBamStat(self, bamfile, chrom, pos, codon_pos):
        bamfile = pysam.AlignmentFile(bamfile, 'rb')
        codon_coverage = 0
        if codon_pos == 0:
            avg_codon_coverage = bamfile.count(chrom, pos-1, pos+2)
        elif codon_pos == 1:
            avg_codon_coverage = bamfile.count(chrom, pos-2, pos+1)
        else:
            avg_codon_coverage = bamfile.count(chrom, pos-3, pos)
        return(avg_codon_coverage)

    def getVarOfInt(self, vcf_path):
        sample_dirs = glob.glob('{0}/*'.format(self.out_path))
        codons_of_int = self.createTable()
        codon_annot = pd.DataFrame()
        for samples in sample_dirs:
            sample_folder = os.path.basename(samples)
            sample_bam = '{0}/output_sorted_RG.bam'.format(samples)
            if not os.path.isdir(samples):
                continue
            sample_vcf = vcf.Reader(filename='{0}/{1}_variants_merged_annotated.vcf'.format(samples, sample_folder))
            sample_name = sample_vcf.samples[0]
            sample_result = []
            sample_annot = []
            sample_vars = {'{0}_{1}_{2}_{3}'.format(var.CHROM, var.POS, var.INFO['RefAA'][0], var.INFO['AltAA'][0]): var.INFO['DP'] for var in sample_vcf}
            for val in codons_of_int.iterrows():
                if '{0}_{1}_{2}_{3}'.format(val[0][0], val[1]['NucPos'], val[1]['Ref'],val[1]['Alt']) in sample_vars.keys():
                    depth = sample_vars['{0}_{1}_{2}_{3}'.format(val[0][0], val[1]['NucPos'],
                                                                           val[1]['Ref'],val[1]['Alt'])]+1
                    sample_result.append(depth)
                    sample_annot.append('{0}->{1}'.format(val[1]['AA'], val[1]['Alt']))
                else:
                    depth = -1 * (self.getBamStat(sample_bam, val[0][0], val[1]['NucPos'],
                                                           val[1]['CodonPos'])+1)
                    sample_result.append(depth)
                    sample_annot.append('WT')
            sample_series = pd.Series(sample_result, index=codons_of_int.index)
            annot_series = pd.Series(sample_annot, index=codons_of_int.index)
            codons_of_int[sample_name] = sample_series
            codon_annot[sample_name] = sample_series

        codons_of_int = codons_of_int.iloc[:,9:].groupby(codons_of_int.index).sum()
        codon_annot = codon_annot.as_matrix()
        sns.set()
        sns.set_context("notebook")
        plt.figure(figsize=(24,30))
        cmap = sns.color_palette("RdBu", n_colors=7)
        heatmap = sns.heatmap(codons_of_int, linewidths=.5,  square=False)
        fig = heatmap.get_figure()
        fig.savefig('{0}/heatmap.png'.format(self.out_path))



if __name__ == '__main__':
    fasta_path = sys.argv[1]
    bed_path = sys.argv[2]
    voi_path = sys.argv[3]
    out_path = sys.argv[4]
    summarizer = Summary(fasta_path, bed_path, voi_path, out_path)
    summarizer.getVarOfInt(out_path)
