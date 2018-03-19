import os
import re
import sys
import glob
import warnings
import numpy as np
import pandas as pd
import logging
import pysam
import pathlib
from pyamd.parsers.vcf import Vcf
from pyamd.parsers.bed import Bed
#from pyamd.readers import Bed

class Summary:

    def __init__(self, fasta, bed, voi, out_path):
        self.fasta = fasta
        self.bed = bed
        self.voi = voi
        self.out_path = out_path

    def getVarOfInt(self):
        if pathlib.Path(self.voi).suffix == '.xlsx':
            voi_table = pd.read_excel(self.voi)
        elif pathlib.Path(self.voi).suffix == '.csv':
            voi_table = pd.read_table(self.voi, sep=',')
        elif pathlib.Path(self.voi).suffix == '.tsv':
            voi_table = pd.read_table(self.voi, sep='\t')
        voi_table['Variant'] = voi_table['Gene'] + ':' + voi_table['RefAA'] + voi_table['AAPos'].astype(str) + voi_table['AltAA']
        voi_table['SNP'] = voi_table['RefAA']+voi_table['AAPos'].astype(str) + voi_table['AltAA']
        voi_table.set_index(['Variant'], inplace=True)
        return(voi_table)

    def getGeneStats(self, bam_path):
        bam_file = pysam.AlignmentFile(bam_path, 'rb')
        bed = Bed(self.bed)
        gene_stats = dict()
        for bed_rec in bed.getExonTable():
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

    def getIntronTables(self):
        vcf_files = glob.glob('{0}/*/*_variants_merged_annotated.vcf'.format(self.out_path))
        vcf_df = pd.DataFrame()
        vcf_dict = {'Gene': [], 'Pos': [], 'Qual': [], 'Ref': [], 'Alt': [], 'AAPos': [],
                    'AltCodon': [], 'RefCodon': [], 'RefAA': [], 'AltAA': [], 'DP': [],
                    'AF': [], 'Conf': [], 'Exon': [], 'Chrom': []}
        vcf_var = list()
        vcf_sample = list()
        vcf_gene = list()
        for files in vcf_files:
            vcf = Vcf.Reader(files)
            vcf_file = vcf.read()
            barcode = re.compile('_[ATGC]*-[ATGC]*')
            sample = barcode.split(vcf.samples[0])[0]
            for var in vcf_file:
                #print(var.INFO)
                if var.INFO['Exon'][0] == 'Intron':
                    vcf_dict['Chrom'].append(var.CHROM)
                    vcf_dict['Gene'].append(var.CHROM)
                    vcf_dict['Pos'].append(var.POS)
                    vcf_dict['Qual'].append(var.QUAL)
                    vcf_dict['Ref'].append(var.REF[0])
                    vcf_dict['Alt'].append(str(var.ALT[0]))
                    vcf_dict['Exon'].append('Intron')
                    vcf_dict['AAPos'].append(np.nan)
                    vcf_dict['RefCodon'].append('NA')
                    vcf_dict['AltCodon'].append('NA')
                    vcf_dict['RefAA'].append('NA')
                    vcf_dict['AltAA'].append('NA')
                    vcf_dict['DP'].append(var.INFO['DP'][0])
                    vcf_dict['AF'].append(float(var.INFO['Freq'][0])*100)
                    vcf_dict['Conf'].append(int(var.INFO['Conf'][0]))
                    vcf_gene.append(var.CHROM)
                    vcf_var.append('{0}:{1}{2}{3}'.format(var.CHROM, var.REF[0], var.POS, str(var.ALT[0])))
                    vcf_sample.append(sample)
        vcf_index = [np.array(vcf_sample), np.array(vcf_var)]
        vcf_df = pd.DataFrame(vcf_dict, index=vcf_index)
        vcf_df.index.names = ['Sample', 'Variant']
        return(vcf_df)


    def getVarTables(self):
        voi_table = self.getVarOfInt()
        vcf_files = glob.glob('{0}/*/*_variants_merged_annotated.vcf'.format(self.out_path))
        vcf_df = pd.DataFrame()
        vcf_dict = {'Gene' : [], 'Pos' : [], 'Qual' : [], 'Ref' : [], 'Alt' : [], 'AAPos' : [], 'RefCodon' : [],
                    'AltCodon' : [], 'RefAA' : [], 'AltAA' : [], 'DP' : [], 'AF' : [], 'Conf': [], 'Exon' : [],
                    'Chrom' : []}
        vcf_var = list()
        vcf_sample = list()
        vcf_gene = list()
        var_sample = list()
        voi_df = self.getVarOfInt()
        #print(voi_df)
        for files in vcf_files:
            vcf = Vcf.Reader(files)
            vcf_file = vcf.read()
            barcode = re.compile('_[ATGC]*-[ATGC]*')
            sample = barcode.split(vcf.samples[0])[0]
            count = 0
            for var in vcf_file:
                #print(var.CHROM, var.POS)
                #print(var.Samples)
                #print(var.INFO)
                if var.CHROM == None or var.INFO['RefAA'][0] == None or var.INFO['AAPos'][0] == None or var.INFO['AltAA'][0] == None:
                    continue
                if '{0}:{1}{2}{3}'.format(var.CHROM, var.INFO['RefAA'][0], var.INFO['AAPos'][0], var.INFO['AltAA'][0]) in voi_df.index:
                    count += 1
                if '{4}{0}:{1}{2}{3}'.format(var.CHROM, var.INFO['RefAA'][0], var.INFO['AAPos'][0], var.INFO['AltAA'][0],sample) in var_sample:
                    index = 0
                    for pos in range(len(vcf_var)):
                        if sample == vcf_sample[pos] and '{0}:{1}{2}{3}'.format(var.CHROM, var.INFO['RefAA'][0], var.INFO['AAPos'][0], var.INFO['AltAA'][0]) == vcf_var[pos]:
                            index = pos
                    vcf_dict['Ref'][index] = '{0},{1}'.format(vcf_dict['Ref'][index], var.REF[0])
                    vcf_dict['Alt'][index] = '{0},{1}'.format(vcf_dict['Alt'][index], str(var.ALT[0]))
                else:
                    vcf_dict['Chrom'].append(var.CHROM)
                    vcf_dict['Gene'].append(var.INFO['Gene'][0])
                    vcf_dict['Pos'].append(var.POS)
                    vcf_dict['Qual'].append(var.QUAL)
                    vcf_dict['Ref'].append(var.REF[0])
                    vcf_dict['Alt'].append(str(var.ALT[0]))
                    vcf_dict['Exon'].append(var.INFO['Exon'][0])
                    vcf_dict['AAPos'].append(int(var.INFO['AAPos'][0]))
                    vcf_dict['RefCodon'].append(var.INFO['RefCodon'][0])
                    vcf_dict['AltCodon'].append(var.INFO['AltCodon'][0])
                    vcf_dict['RefAA'].append(var.INFO['RefAA'][0])
                    vcf_dict['AltAA'].append(var.INFO['AltAA'][0])
                    vcf_dict['DP'].append(var.INFO['DP'][0])
                    vcf_dict['AF'].append(float(var.INFO['Freq'][0]) * 100)
                    vcf_dict['Conf'].append(int(var.INFO['Conf'][0]))
                    vcf_gene.append(var.CHROM)
                    vcf_var.append('{0}:{1}{2}{3}'.format(var.INFO['Gene'][0], var.INFO['RefAA'][0], var.INFO['AAPos'][0], var.INFO['AltAA'][0]))
                    vcf_sample.append(sample)
                    if var.INFO['DP'][0] > 0:
                        var_sample.append('{4}{0}:{1}{2}{3}'.format(var.INFO['Gene'][0], var.INFO['RefAA'][0], var.INFO['AAPos'][0], var.INFO['AltAA'][0], sample))
                    else:
                        var_sample.append('{4}{0}:{1}{2}NA'.format(var.INFO['Gene'][0], var.INFO['RefAA'][0], var.INFO['AAPos'][0], var.INFO['AltAA'][0], sample))
                    #count += 1


            if count == 0:
                for variants, rec in voi_df.iterrows():
                    vcf_dict['Chrom'].append(rec.Chrom)
                    vcf_dict['Gene'].append(rec.Gene)
                    vcf_dict['Pos'].append(np.nan)
                    vcf_dict['Qual'].append(np.nan)
                    vcf_dict['Ref'].append(np.nan)
                    vcf_dict['Alt'].append(np.nan)
                    vcf_dict['Exon'].append(np.nan)
                    vcf_dict['AAPos'].append(rec.AAPos)
                    vcf_dict['RefCodon'].append(np.nan)
                    vcf_dict['AltCodon'].append(np.nan)
                    vcf_dict['RefAA'].append(rec.RefAA)
                    vcf_dict['AltAA'].append(rec.AltAA)
                    vcf_dict['DP'].append(0)
                    vcf_dict['AF'].append(np.nan)
                    vcf_dict['Conf'].append(2)
                    #print('{0}{1}{2}'.format(var.group('RefAA'), var.group('AAPos'), var.group('RefAA')))
                    vcf_var.append(variants)
                    vcf_sample.append(sample)

        vcf_index = [np.array(vcf_sample), np.array(vcf_var)]
        vcf_df = pd.DataFrame(vcf_dict, index=vcf_index)
        vcf_df.index.names = ['Sample', 'Variant']
        #print(vcf_df.head())
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
        exp_voi['FinalCall'] = exp_voi['SNP']
        for index, series in exp_voi.iterrows():
            #print(exp_voi.at[index, 'FinalCall'])
            if pd.isnull(series['DP']):
                exp_voi.at[index, 'FinalCall'] = 'WT'
            elif pd.isnull(series['Alt']):
                var_reg = re.match(r'(?P<RefAA>[DTSEPGACVMILYFHKRWQN])(?P<AAPos>\d+)(?P<AltAA>[DTSEPGACVMILYFHKRWQN])', series['SNP'])
                exp_voi.at[index, 'FinalCall'] = '{0}{1}{0}'.format(var_reg.group('RefAA'), var_reg.group('AAPos'))
                #print('{0}{1}{0}'.format(var_reg.group('RefAA'), var_reg.group('AAPos')))
                #print(series['FinalCall'])
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
        return(exp_nov)

    def getBamStat(self, bamfile, chrom, start, stop):
        bamfile = pysam.AlignmentFile(bamfile, 'rb')
        avg_codon_coverage = bamfile.count(chrom, start, stop)
        return(avg_codon_coverage)

    def getNucPos(self, gene, aapos):
        bed = Bed(self.bed)
        exon_table = bed.getExonTable()
        bed_list = list()
        for records in exon_table:
            if gene == records.gene:
                bed_list += [val for val in range(records.start, records.stop+1)]
        bed_list = [bed_list[ind:ind+3] for ind in range(0, len(bed_list),3)]
        try:
            return(bed_list[int(aapos)-1])
        except ValueError:
            return(np.nan)

    def getDepthStats(self, var_df):
        depth_list = list()
        for row, value in var_df.iterrows():
            bamfile = glob.glob('{0}/{1}*/output_fixmate_sorted_RG.bam'.format(self.out_path, row[0]))[0]
            nuc_pos = self.getNucPos(value.Gene_y, value.AAPos_y)
            if nuc_pos == np.nan:
                nuc_pos = [value.Pos -1, value.Pos + 1]
            depth = self.getBamStat(bamfile, value.Chrom_y, nuc_pos[0], nuc_pos[1])
            depth_list.append(depth)            #np.log10(depth+1))
        var_df['DP'] = pd.Series(depth_list, index=var_df.index)
        return(var_df)

    def getNovDepthStats(self, var_df):
        depth_list = list()
        for row, value in var_df.iterrows():
            bamfile = glob.glob('{0}/{1}*/output_fixmate_sorted_RG.bam'.format(self.out_path, row[0]))[0]
            nuc_pos = self.getNucPos(value.Gene, value.AAPos)
            if nuc_pos == np.nan:
                nuc_pos = [value.Pos -1, value.Pos + 1]
            depth = self.getBamStat(bamfile, value.Chrom, nuc_pos[0], nuc_pos[1])
            depth_list.append(depth)
        var_df['DP'] = pd.Series(depth_list, index=var_df.index)
        return(var_df)

    def getVarStats(self, vcf_file):
        vcf_file = Vcf.Reader(vcf_file)
        vcf_reader = vcf_file.read()
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
        for variant in vcf_reader:
            total += 1
            if variant.INFO['Exon'][0] == 'Intron':
                intronic += 1
            else:
                exonic += 1
                if variant.INFO['Conf'][0] == 2:
                    verfied += 1
                if variant.INFO['RefAA'][0] == variant.INFO['AltAA'][0]:
                    syn += 1
                else:
                    nsyn += 1
                if '{0}{1}'.format(variant.REF[0], str(variant.ALT[0])) in trasition:
                    trans += 1
                else:
                    tranv += 1
        return(total, verfied, exonic, intronic, syn, nsyn, trans, tranv)


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
