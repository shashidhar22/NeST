import os
import re
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

    def plotHeatMap(self, codon_of_int):
        sns.set()
        sns.set_context('notebook')
        plt.figure(figsize=(24,30))
        cmap = sns.diverging_palette(0,359,sep=90, as_cmap=True)
        heatmap = sns.heatmap(codon_of_int, linewidths=0.5, cmap="Blues") #sns.cubehelix_palette(8, start=.5, rot=-.75, as_cmap=True))
        fig = heatmap.get_figure()
        fig.savefig('{0}/heatmap.png'.format(self.out_path))
        return


    def getVarOfInt(self, vcf_path):
        sample_dirs = glob.glob('{0}/*'.format(self.out_path))
        codons_of_int = self.createTable()
        codons_of_int.to_excel('{0}/Codons_of_interest.xlsx'.format(self.out_path))
        for samples in sample_dirs:
            sample_folder = os.path.basename(samples)
            sample_bam = '{0}/output_sorted_RG.bam'.format(samples)
            if not os.path.isdir(samples):
                continue
            sample_vcf = vcf.Reader(filename='{0}/{1}_variants_merged_annotated.vcf'.format(samples, sample_folder))
            sample_name = sample_vcf.samples[0]
            sample_result = []
            sample_annot = []
            sample_vars = {'{0}_{1}_{2}_{3}'.format(var.CHROM, var.POS, var.INFO['RefAA'][0], var.INFO['AltAA'][0]): float(var.INFO['AlFreq'][0]) * 100.0 for var in sample_vcf}
            for val in codons_of_int.iterrows():
                if '{0}_{1}_{2}_{3}'.format(val[0][0], val[1]['NucPos'], val[1]['Ref'],val[1]['Alt']) in sample_vars.keys():
                    depth = sample_vars['{0}_{1}_{2}_{3}'.format(val[0][0], val[1]['NucPos'],

                                                                           val[1]['Ref'],val[1]['Alt'])]
                else:
                    depth = 0.0 

                sample_result.append(depth)
                            
            sample_series = pd.Series(sample_result, index=codons_of_int.index)
            codons_of_int[sample_name] = sample_series

        codons_of_int.to_excel('{0}/Codons_of_interest.xlsx'.format(self.out_path))
        codons_of_int = codons_of_int.iloc[:,9:].groupby(codons_of_int.index).max()
        codons_of_int['label'] = pd.Series([var[0] + ':' + str(var[1]) for var in codons_of_int.index.tolist()], index=codons_of_int.index)
        codons_of_int.set_index('label', inplace=True)
        #del codons_of_int['Gene']
        codons_of_int.to_excel('{0}/al_freq_table.xlsx'.format(self.out_path))
        self.plotHeatMap(codons_of_int)
        return

    def getAfHeatmap(self, vcf_path):
        voi_table = pd.read_excel(self.voi, sheetname=1)
        experiment_df = pd.DataFrame()
        for vcf_files in glob.glob('{0}/*/*_variants_merged_annotated.vcf'.format(self.out_path)):
            vcf_dict = {'Gene' : [], 'Pos' : [], 'Qual' : [], 'Ref' : [], 'Alt' : [], 'CodonPos' : [], 'RefCodon' : [], 'AltCodon' : [], 'RefAA' : [], 'AltAA' : [], 'DP' : [], 'AF' : [], 'Conf': [], 'Exon' : []}
            vcf_index = list()
            vcf_file = vcf.Reader(filename=vcf_files)
            barcode = re.compile('_[ATGC]*-[ATGC]*')
            sample = barcode.split(vcf_file.samples[0])[0]
            for variant in vcf_file:
                vcf_dict['Gene'].append(variant.CHROM)
                vcf_dict['Pos'].append(variant.POS)
                vcf_dict['Qual'].append(variant.QUAL)
                vcf_dict['Ref'].append(variant.REF)
                vcf_dict['Alt'].append(str(variant.ALT[0]))
                vcf_dict['Exon'].append(variant.INFO['ExonNumber'][0])
                vcf_dict['CodonPos'].append(np.nan if variant.INFO['CodonPos'][0] == 'NA' else int(variant.INFO['CodonPos'][0]))
                vcf_dict['RefCodon'].append(np.nan if variant.INFO['RefCodon'][0] == 'NA' else variant.INFO['RefCodon'][0])
                vcf_dict['AltCodon'].append(np.nan if variant.INFO['AltCodon'][0] == 'NA' else variant.INFO['AltCodon'][0])
                vcf_dict['RefAA'].append(np.nan if variant.INFO['RefAA'][0] == 'NA' else variant.INFO['RefAA'][0])
                vcf_dict['AltAA'].append(np.nan if variant.INFO['AltAA'][0] == 'NA' else variant.INFO['AltAA'][0])
                vcf_dict['DP'].append(variant.INFO['DP'])
                vcf_dict['AF'].append(float(variant.INFO['AlFreq'][0]) * 100)
                vcf_dict['Conf'].append(int(variant.INFO['Found'][0]))
                vcf_index.append(sample)
            vcf_df = pd.DataFrame(vcf_dict, index=vcf_index)
            vcf_exon = vcf_df[vcf_df['Exon'] != 'Intron']
            vcf_exon['SNP'] = vcf_exon['RefAA'] + vcf_exon['CodonPos'].map(int).map(str) + vcf_exon['AltAA']
            voi_exon = vcf_exon.merge(voi_table, on=['Gene', 'SNP'], how='right')
            sample = [vcf_index[0]] *  len(voi_exon)
            voi_exon['Sample'] = pd.Series(sample, index=voi_exon.index)
            voi_exon.set_index('Sample', inplace=True)
            experiment_df = experiment_df.append(voi_exon)
        experiment_df.replace({'PfCRT':'CRT', 'PfMDR1':'MDR1'}, inplace=True)
        experiment_df['SNP'] = experiment_df['Gene'] + ':' + experiment_df['SNP'] 
        experiment_af = experiment_df.pivot(experiment_df.index, 'SNP')['AF'].transpose()
        af_mask = experiment_af.isnull()
        sns.set()
        sns.set_style('whitegrid')
        fig, ax = plt.subplots()
        fig.set_size_inches(24, 24)
        cbar_ax = fig.add_axes([.92, .3, .02, .4])
        heatmap_af = sns.heatmap(experiment_af, linewidths=0.5, vmin=0.0, cmap="Blues", ax=ax, cbar_ax=cbar_ax, mask=af_mask, square=True)
        fig_af = heatmap_af.get_figure()
        fig_af.savefig('{0}/alfreq_heatmap.png'.format(self.out_path))
        experiment_count = experiment_af.count(axis=1, numeric_only=True).to_frame('Count')
        experiment_count['Var'] = experiment_count.index
        experiment_count['Gene'], experiment_count['SNP'] = experiment_count['Var'].str.split(':',1).str
        del experiment_count['Var']
        sns.set(font_scale=2)
        plt.figure(figsize=(20, 20))
        stripplot = sns.stripplot(y=experiment_af.index, x=experiment_af.count(axis=1, numeric_only=True), size=15, color='black')
        plots = stripplot.get_figure()
        plots.savefig('{0}/frequency.png'.format(self.out_path))
		#Plot depth heatmap
        bed_file = pd.read_table(self.bed, header=None, sep='\t', names=['Gene','Start','Stop','Exon','Info','Strand'] )
        bed_file.replace({'PfCRT':'CRT', 'PfMDR1':'MDR1'}, inplace=True)
        cdna_dict = dict()
        for index, row in bed_file.iterrows():
            for base in range(row['Start'], row['Stop']):
                try:
                    cdna_dict[row['Gene']].append(base)
                except KeyError:
                    cdna_dict[row['Gene']] = [base]

        for key, values in cdna_dict.items():
            values = [values[i:i+3] for i in range(0, len(values),3)]
            cdna_dict[key] = values
        pfs = {'CRT' : 'PfCRT', 'MDR1' : 'PfMDR1'}
        experiment_df['VOIPos'] = experiment_df['SNP'].str.extract('(\d+)').astype(int)
        for index, row in experiment_df.iterrows():
            filename = glob.glob('{0}/{1}_*/output_sorted_RG.bam'.format(self.out_path, index))[0]
            sample_bam = pysam.AlignmentFile(filename, 'rb')
            if pd.isnull(row['DP']):
                try:
                    avg_codon_coverage = sample_bam.count(pfs[row['Gene']], min(cdna_dict[row['Gene']][row['VOIPos']]), max(cdna_dict[row['Gene']][row['VOIPos']]))
                except KeyError:
                    avg_codon_coverage = sample_bam.count(row['Gene'], min(cdna_dict[row['Gene']][row['VOIPos']]), max(cdna_dict[row['Gene']][row['VOIPos']]))
                experiment_df['DP'].loc[index] = avg_codon_coverage
        experiment_dp = experiment_df.pivot(experiment_df.index, 'SNP')['DP'].transpose()
        #experiment_dp = np.log(experiment_dp + 1)
        dp_mask = experiment_dp.isnull()
        sns.set()
        sns.set_style('whitegrid')
        #sns.set_context('notebook')
        fig, ax = plt.subplots()
        fig.set_size_inches(24, 24)
        cbar_ax = fig.add_axes([.92, .3, .02, .4])
        heatmap_dp = sns.heatmap(experiment_dp, linewidths=0.5, vmin=0.0, cmap="Blues", ax=ax, cbar_ax=cbar_ax, mask=dp_mask, square=True)
        fig_dp = heatmap_dp.get_figure()
        fig_dp.savefig('{0}/depth_heatmap.png'.format(self.out_path))
        return

if __name__ == '__main__':
    fasta_path = sys.argv[1]
    bed_path = sys.argv[2]
    voi_path = sys.argv[3]
    out_path = sys.argv[4]
    summarizer = Summary(fasta_path, bed_path, voi_path, out_path)
    summarizer.getVarOfInt(out_path)
