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
import warnings
import seaborn as sns

warnings.filterwarnings('ignore')

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


    def vcfToTable(self, vcf_path):
        voi_table = pd.read_excel(self.voi, sheetname=1)
        experiment_df = pd.DataFrame()
        experiment_nov_df = pd.DataFrame()
        for vcf_files in glob.glob('{0}/*/*_variants_merged_annotated.vcf'.format(self.out_path)):
            vcf_dict = {'Gene' : [], 'Pos' : [], 'Qual' : [], 'Ref' : [], 'Alt' : [], 'CodonPos' : [], 'RefCodon' : [],
                        'AltCodon' : [], 'RefAA' : [], 'AltAA' : [], 'DP' : [], 'AF' : [], 'Conf': [], 'Exon' : []}
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
            vcf_exon.to_excel('{0}/vcf_exon_table.xlsx'.format(self.out_path))
            voi_table.to_excel('{0}/voi_var_table.xlsx'.format(self.out_path))
            #nov_exon = vcf_exon[(vcf_exon['SNP'] != voi_table['SNP']) & (vcf_exon['Gene'] == voi_table['Gene'])]
            voi_exon = vcf_exon.merge(voi_table, on=['Gene', 'SNP'], how='right')
            #nov_exon = vcf_exon[(vcf_exon['SNP'] != voi_table['SNP']) & (vcf_exon['Gene'] == voi_table['Gene'])]
            sample = [vcf_index[0]] *  len(voi_exon)
            voi_exon['Sample'] = pd.Series(sample, index=voi_exon.index)
            voi_exon.set_index('Sample', inplace=True)
            experiment_df = experiment_df.append(voi_exon)
            #experiment_nov_df = experiment_nov_df.append(nov_exon)
        return(experiment_df)

    def getVarTable(self, vcf_path):
        #Iterate through all vcf files and generate table for variants of interest and novel variants
        experiment_df = self.vcfToTable(vcf_path)
        experiment_df.replace({'PfCRT':'CRT', 'PfMDR1':'MDR1'}, inplace=True)
        #experiment_nov_df.replace({'PfCRT':'CRT', 'PfMDR1':'MDR1'}, inplace=True)
        experiment_df['SNP'] = experiment_df['Gene'] + ':' + experiment_df['SNP']
        #experiment_nov_df['SNP'] = experiment_nov_df['Gene'] + ':' + experiment_nov_df['SNP']
        #Create allele frequency table for all samples
        experiment_af = experiment_df.pivot(experiment_df.index, 'SNP')['AF'].transpose()
        #experiment_nov_af = experiment_nov_df.pivot(experiment_nov_df.index, 'SNP')['AF'].transpose()
        #Create count table containing number of samples containing a particular SNP
        experiment_count = experiment_af.count(axis=1, numeric_only=True).to_frame('Count')
        #experiment_nov_count = experiment_nov_af.count(axis=1, numeric_only=True).to_frame('Count')
        experiment_count['Var'] = experiment_count.index
        experiment_count['Gene'], experiment_count['SNP'] = experiment_count['Var'].str.split(':',1).str
        del experiment_count['Var']
        #Create depth of coverage table for all samples
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
        #experiment_nov_df['VOIPos'] = experiment_nov_df['SNP'].str.extract('(\d+)').astype(int)
        #Get codon depth for vairiants of interest
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
        #Get codon depth for novel vairiants
        #for index, row in experiment_nov_df.iterrows():
        #    filename = glob.glob('{0}/{1}_*/output_sorted_RG.bam'.format(self.out_path, index))[0]
        #    sample_bam = pysam.AlignmentFile(filename, 'rb')
        #    if pd.isnull(row['DP']):
        #        try:
        #            avg_codon_coverage = sample_bam.count(pfs[row['Gene']], min(cdna_dict[row['Gene']][row['VOIPos']]), max(cdna_dict[row['Gene']][row['VOIPos']]))
        #        except KeyError:
        #            avg_codon_coverage = sample_bam.count(row['Gene'], min(cdna_dict[row['Gene']][row['VOIPos']]), max(cdna_dict[row['Gene']][row['VOIPos']]))
        #        experiment_nov_df['DP'].loc[index] = avg_codon_coverage
        #experiment_nov_dp = experiment_nov_df.pivot(experiment_df.index, 'SNP')['DP'].transpose()
        return(experiment_df, experiment_af, experiment_count, experiment_dp)

    def plotHeatMap(self, data_frame, title, mask):
        dp_mask = mask
        sns.set()
        sns.set_style('whitegrid')
        fig, ax = plt.subplots()
        fig.set_size_inches(24, 24)
        cbar_ax = fig.add_axes([.92, .3, .02, .4])
        heatmap_dp = sns.heatmap(data_frame, linewidths=0.5, vmin=0.0, cmap="Blues", ax=ax, cbar_ax=cbar_ax, mask=dp_mask, square=True)
        fig_dp = heatmap_dp.get_figure()
        fig_dp.savefig('{0}/{1}_heatmap.png'.format(self.out_path, title))
        return

    def plotCountPlot(self, data_frame, title):
        sns.set(font_scale=2)
        plt.figure(figsize=(20, 20))
        stripplot = sns.stripplot(y=data_frame.index, x=data_frame.count(axis=1, numeric_only=True), size=15, color='black')
        plots = stripplot.get_figure()
        plots.savefig('{0}/{1}_frequency.png'.format(self.out_path, title))


    def getSummary(self, voi_df, voi_af, voi_count, voi_dp):
        #Create masks for heatmap
        dp_voi_mask = voi_dp.isnull()
        af_voi_mask = voi_af.isnull()
        #dp_nov_mask = nov_dp.isnull()
        #af_nov_mask = nov_af.isnull()
        #Plot depth heatmap for variants of interest
        self.plotHeatMap(voi_dp, 'voi_depth', dp_voi_mask)
        #Plot depth heatmap for novel variants
        #self.plotHeatMap(nov_dp, 'nov_depth', dp_nov_mask)
        #Plot allele frequency heatmap for variants of interest
        self.plotHeatMap(voi_dp, 'voi_alfreq', af_voi_mask)
        #plot allele frequency heatmap for novel variants
        #self.plotHeatMap(nov_dp, 'nov_alfreq', af_voi_mask)
        #Plot frequency plots
        self.plotCountPlot(voi_af, 'voi')
        #self.plotCountPlot(nov_af, 'nov')
        return

if __name__ == '__main__':
    fasta_path = sys.argv[1]
    bed_path = sys.argv[2]
    voi_path = sys.argv[3]
    out_path = sys.argv[4]
    summarizer = Summary(fasta_path, bed_path, voi_path, out_path)
    summarizer.getVarOfInt(out_path)
