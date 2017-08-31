import os
import re
import sys
import vcf
import glob
import math
import pandas as pd
import pysam
import logging
from pyamd.parsers import Fasta
from pyamd.readers import Bed
from pyamd.annotater import Annotate
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import warnings
import seaborn as sns

print('matplotlib', matplotlib.__version__)
print('numpy', np.__version__)
print('pysam', pysam.__version__)
print('seaborn', sns.__version__)
print('pandas', pd.__version__)
warnings.filterwarnings('ignore')

class Summary:

    def __init__(self, fasta, bed, voi, out_path):
        self.fasta = fasta
        self.bed = bed
        self.voi = voi
        self.out_path = out_path
        self.logger = logging.getLogger('Mars.sample_runner.summarize')

    def getVarOfInt(self):
        voi_table = pd.read_excel(self.voi, sheetname=1, parse_cols="B:C")
        voi_df = voi_table['SNP'].str.extract('(?P<RefAA>[a-zA-Z]?)(?P<AAPos>[0-9]*)(?P<AltAA>[a-zA-Z]?)', expand=True)
        voi_df['Gene'] = voi_table['Gene']
        voi_df['SNP'] = voi_table['SNP']
        voi_df['Variant'] = voi_df['Gene'] + ':' + voi_df['SNP']
        voi_df.set_index(['Variant'], inplace=True)
        return(voi_df)


    def getGeneStats(self, bam_path):
        bam_file = pysam.AlignmentFile(bam_path, 'rb')
        bed = Bed(self.bed)
        gene_stats = dict()
        for bed_rec in bed.read():
            try:
                gene_stats[bed_rec.chrom] += list(np.sum(bam_file.count_coverage(bed_rec.chrom, bed_rec.start, bed_rec.stop+1), axis=0))
            except KeyError:
                gene_stats[bed_rec.chrom] = list(np.sum(bam_file.count_coverage(bed_rec.chrom, bed_rec.start, bed_rec.stop+1), axis=0))
        return(gene_stats)

    def getExpCoverage(self):
        bam_files = glob.glob('{0}/*/output_sorted_RG.bam'.format(self.out_path))
        barcode = re.compile('_[ATGC]*-[ATGC]*')
        sample_gstat = dict()
        for files in bam_files:
            sample_dir = os.path.basename(os.path.dirname(files))
            sample = barcode.split(sample_dir)[0]
            sample_gstat[sample] = self.getGeneStats(files)
        return(sample_gstat)

    def checkDepthPass(self):
        exp_depth = self.getExpCoverage()
        exp_pass = dict()
        for samples in exp_depth:
            gene_pass = dict()
            for genes in exp_depth[samples]:
                if np.percentile(exp_depth[samples][genes], 50) < 10:
                    gene_pass[genes] = False
                else:
                    gene_pass[genes] = True
            exp_pass[samples] = gene_pass
        return(exp_pass)

    def getVarTables(self):
        voi_table = self.getVarOfInt()
        vcf_files = glob.glob('{0}/*/*_variants_merged_annotated.vcf'.format(self.out_path))
        vcf_df = pd.DataFrame()
        vcf_dict = {'Gene' : [], 'Pos' : [], 'Qual' : [], 'Ref' : [], 'Alt' : [], 'CodonPos' : [], 'RefCodon' : [],
                    'AltCodon' : [], 'RefAA' : [], 'AltAA' : [], 'DP' : [], 'AF' : [], 'Conf': [], 'Exon' : []}
        vcf_var = list()
        vcf_sample = list()
        vcf_gene = list()
        voi_df = self.getVarOfInt()
        for files in vcf_files:
            vcf_file = vcf.Reader(filename=files)
            barcode = re.compile('_[ATGC]*-[ATGC]*')
            sample = barcode.split(vcf_file.samples[0])[0]
            count = 0
            for var in vcf_file:
                if var.CHROM == 'NA': #var.INFO['RefAA'][0] == 'NA' or var.INFO['CodonPos'][0] == 'NA' or var.INFO['AltAA'][0] == 'NA':
                    continue
                if '{0}:{1}{2}{3}'.format(var.CHROM, var.INFO['RefAA'][0], var.INFO['CodonPos'][0], var.INFO['AltAA'][0]) in voi_df.index:
                    count += 1
                vcf_dict['Gene'].append(var.CHROM)
                vcf_dict['Pos'].append(var.POS)
                vcf_dict['Qual'].append(var.QUAL)
                vcf_dict['Ref'].append(var.REF)
                vcf_dict['Alt'].append(str(var.ALT[0]))
                vcf_dict['Exon'].append(var.INFO['ExonNumber'][0])
                vcf_dict['CodonPos'].append(int(var.INFO['CodonPos'][0]))
                vcf_dict['RefCodon'].append(var.INFO['RefCodon'][0])
                vcf_dict['AltCodon'].append(var.INFO['AltCodon'][0])
                vcf_dict['RefAA'].append(var.INFO['RefAA'][0])
                vcf_dict['AltAA'].append(var.INFO['AltAA'][0])
                vcf_dict['DP'].append(var.INFO['DP'])
                vcf_dict['AF'].append(float(var.INFO['AlFreq'][0]) * 100)
                vcf_dict['Conf'].append(int(var.INFO['Found'][0]))
                vcf_gene.append(var.CHROM)
                vcf_var.append('{0}:{1}{2}{3}'.format(var.CHROM, var.INFO['RefAA'][0], var.INFO['CodonPos'][0], var.INFO['AltAA'][0]))
                vcf_sample.append(sample)
                #count += 1
            if count == 0:
                self.logger.info('No variants found; adding ref calls to dataframe')
                for variants, rec in voi_df.iterrows():
                    vcf_dict['Gene'].append(rec.Gene)
                    vcf_dict['Pos'].append(np.nan)
                    vcf_dict['Qual'].append(np.nan)
                    vcf_dict['Ref'].append(np.nan)
                    vcf_dict['Alt'].append(np.nan)
                    vcf_dict['Exon'].append(np.nan)
                    vcf_dict['CodonPos'].append(rec.AAPos)
                    vcf_dict['RefCodon'].append(np.nan)
                    vcf_dict['AltCodon'].append(np.nan)
                    vcf_dict['RefAA'].append(rec.RefAA)
                    vcf_dict['AltAA'].append(rec.AltAA)
                    vcf_dict['DP'].append(0)
                    vcf_dict['AF'].append(np.nan)
                    vcf_dict['Conf'].append(2)
                    vcf_var.append(variants)
                    vcf_sample.append(sample)
        vcf_index = [np.array(vcf_sample), np.array(vcf_var)]
        vcf_df = pd.DataFrame(vcf_dict, index=vcf_index)
        vcf_df.index.names = ['Sample', 'Variant']
        vcf_df.to_excel('Varcalls_table.xlsx')
        return(vcf_df)

    def getRepSnps(self):
        exp_df = self.getVarTables()
        voi_df = self.getVarOfInt()
        samp = 0
        exp_voi = pd.DataFrame()
        for sample, var_df in exp_df.groupby(level=0):
            sam_index = list()
            var_df = var_df.reset_index(level=0)
            var_voi = var_df.merge(voi_df, how='right', left_index=True, right_index=True)
            sam_index = [sample] * len(var_voi)
            var_index = [np.array(sam_index), np.array(var_voi.index)]
            var_voi.set_index(var_index, inplace=True)
            var_voi.index.names = ['Sample', 'Variant']
            exp_voi = exp_voi.append(var_voi)
            #print(var_voi.head())
        return(exp_voi)

    def getNovSnps(self):
        exp_df = self.getVarTables()
        voi_df = self.getVarOfInt()
        exp_nov = pd.DataFrame()
        for sample, var_df  in exp_df.groupby(level=0):
            sam_index = list()
            var_df = var_df.reset_index(level=0)
            var_nov = var_df[~var_df.index.isin(voi_df.index)]
            sam_index = [sample] * len(var_nov)
            var_index = [np.array(sam_index), np.array(var_nov.index)]
            var_nov.set_index(var_index, inplace=True)
            var_nov.index.names = ['Sample', 'Variant']
            exp_nov = exp_nov.append(var_nov)
        exp_nov = exp_nov[exp_nov.Conf == 2]
        exp_nov.to_excel('novel_variants.xlsx')
        return(exp_nov)

    def getBamStat(self, bamfile, chrom, start, stop):
        bamfile = pysam.AlignmentFile(bamfile, 'rb')
        codon_coverage = 0
        if codon_pos == 0:
            avg_codon_coverage = bamfile.count(chrom, pos-1, pos+2)
        elif codon_pos == 1:
            avg_codon_coverage = bamfile.count(chrom, pos-2, pos+1)
        else:
            avg_codon_coverage = bamfile.count(chrom, pos-3, pos)
        return(avg_codon_coverage)

    def getNucPos(self, gene, aapos):
        bed = Bed(self.bed)
        bed_list = list()
        for records in bed.read():
            if gene == records.chrom:
                bed_list += [val for val in range(records.start, records.stop+1)]
        bed_list = [bed_list[ind:ind+3] for ind in range(0, len(bed_list),3)]
        return(bed_list[int(aapos)-1])

    def getDepthStats(self, var_df):
        depth_list = list()
        for row, value in var_df.iterrows():
            bamfile = glob.glob('{0}/{1}*/output_sorted_RG.bam'.format(self.out_path, row[0]))[0]
            nuc_pos = self.getNucPos(value.Gene_y, value.AAPos)
            depth = self.getBamStat(bamfile, value.Gene_y, nuc_pos[0], nuc_pos[1])
            depth_list.append(np.log10(depth+1))
        var_df['DP'] = pd.Series(depth_list, index=var_df.index)
        return(var_df)

    def getNovDepthStats(self, var_df):
        depth_list = list()
        for row, value in var_df.iterrows():
            bamfile = glob.glob('{0}/{1}*/output_sorted_RG.bam'.format(self.out_path, row[0]))[0]
            nuc_pos = self.getNucPos(value.Gene, value.CodonPos)
            depth = self.getBamStat(bamfile, value.Gene, nuc_pos[0], nuc_pos[1])
            depth_list.append(depth)
        var_df['DP'] = pd.Series(depth_list, index=var_df.index)
        return(var_df)

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


    def plotHeatMap(self, data_frame, title, mask):
        sns.set()
        sns.set_style('whitegrid')
        fig, ax = plt.subplots()
        fig.set_size_inches(24, 24)
        cbar_ax = fig.add_axes([.92, .3, .02, .4])
        if 'af' in title:
            heatmap_dp = sns.heatmap(data_frame, linewidths=0.5, vmin=0.0, vmax=100.0,
                                    cmap="Blues", ax=ax, cbar_ax=cbar_ax,
                                    mask=mask, square=True, linecolor="black")
        else:
            heatmap_dp = sns.heatmap(data_frame, linewidths=0.5, vmin=0.0,
                                    cmap="Blues", ax=ax, cbar_ax=cbar_ax,
                                    mask=mask, square=True, linecolor="black")
        fig_dp = heatmap_dp.get_figure()
        fig_dp.savefig('{0}/{1}_heatmap.png'.format(self.out_path, title))
        return

    def plotCountPlot(self, data_frame, title):
        sns.set(font_scale=2)
        sns.set_style('whitegrid')
        plt.figure(figsize=(20, 20))
        stripplot = sns.stripplot(y=data_frame.index, x=data_frame.count(axis=1, numeric_only=float), size=15, color='black')
        plots = stripplot.get_figure()
        plots.savefig('{0}/{1}_frequency.png'.format(self.out_path, title))


    def getHeatmap(self, voi_df, voi_af, voi_count, voi_dp, nov_df, nov_af, nov_count, nov_dp):
        #Create masks for heatmap
        dp_voi_mask = voi_dp.isnull()
        af_voi_mask = voi_af.isnull()
        self.plotHeatMap(voi_dp, 'voi_depth', dp_voi_mask)
        self.plotHeatMap(voi_dp, 'voi_alfreq', af_voi_mask)
        self.plotCountPlot(voi_af, 'voi')
        return

if __name__ == '__main__':
    fasta_path = sys.argv[1]
    bed_path = sys.argv[2]
    voi_path = sys.argv[3]
    out_path = sys.argv[4]
    summarizer = Summary(fasta_path, bed_path, voi_path, out_path)
    exp_voi = summarizer.getRepSnps()
    #exp_voi.to_excel('variants_of_interest.xlsx')
    exp_voi = summarizer.getDepthStats(exp_voi)
    exp_voi = exp_voi.reset_index(level=1)
    exp_af =  exp_voi.pivot(exp_voi.index, 'Variant')['AF'].transpose()
    exp_af.to_excel('variants_of_interest_af.xlsx')
    exp_dp =  exp_voi.pivot(exp_voi.index, 'Variant')['DP'].transpose()
    exp_dp.to_excel('variants_of_interest_dp2.xlsx')
    exp_nov = summarizer.getNovSnps()
    #exp_nov = summarizer.getDepthStats(exp_nov)
    exp_nov = exp_nov.reset_index(level=1)
    exp_nov_af = exp_nov.pivot(exp_nov.index, 'Variant')['AF'].transpose()
    exp_nov_af.to_excel('novel_variants_af.xlsx')
    exp_nov_af_mask = exp_nov_af.isnull()
    summarizer.plotHeatMap(exp_nov_af, 'nov_af', exp_nov_af_mask)
    depth_pass = summarizer.checkDepthPass()
    for samples in depth_pass:
        print(samples, depth_pass[samples])
