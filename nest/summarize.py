import os
import re
import sys
import gzip
import glob
import json
import pysam
import pathlib
import logging
import warnings
import subprocess
import numpy as np
import pandas as pd
from datetime import datetime
from nest.parsers.vcfReader import Reader
from nest.parsers.bed import Bed
from collections import namedtuple
from collections import OrderedDict

class Summary:
    """ The summary class is built as a part of the NeST framework. The class
    implements modules that can be used to generate summary tables and figures
    for a given study. The module iterates through VCF files of all samples and
    creates tables broken into three categories, Known, Novel Exonic and Novel
    Intronic. It also generates a allele frequency and depth of coverage view
    for the study, which are summarized in table and figuresself."""

    def __init__(self, fasta, bed, voi, out_path):
        """
        Key attributes:
        fasta (str) - Path to reference Fasta file
        bed (str) - Path to reference BED file
        voi (str) - Path to variant of interest file
        out_path (str) - Path to study output directory
        config (dict) - Study config dictionary, this can be obtained by calling
        prepInputs on the input directory
        """
        self.fasta = fasta
        self.bed = bed
        self.voi = voi
        self.out_path = out_path
        self.summary_path = os.path.dirname(__file__)
        self.nest_path = os.path.dirname(self.summary_path)
        self.logger = logging.getLogger('NeST.Summary')

    def getBaseRange(self, chrom, gene, pos):
        """Given the chromosome, gene and amino acid position, it returns the
        genomic range of the codons."""
        bed_reader = Bed(self.bed).getExonTable()
        for region in bed_reader:
            if region.chrom == chrom:
                if region.gene == gene:
                    cds_pos = pos * 3
                    if cds_pos in range(region.exonStart, region.exonStop+1):
                        if region.strand == '+':
                            start = region.start + region.overHang
                            stop = region.stop
                            codon = 1 + region.aaCount
                            for base in range(start, stop, 3):
                                if codon == pos:
                                    return(range(base , base +3 ))
                                else:
                                    codon += 1
                        else:
                            start = region.stop - region.overHang
                            stop = region.start
                            codon = 1 + region.aaCount
                            for base in range(start, stop, -3):
                                if codon == pos:
                                    return(range(base -2, base +1))
                                else:
                                    codon += 1

    def getVarOfInt(self):
        """Return dataframe of variants of interest."""
        if self.voi == None:
            return(None)
        if os.path.splitext(self.voi)[1] == '.xlsx':
            voi_table = pd.read_excel(self.voi)
        elif os.path.splitext(self.voi)[1] == '.csv':
            voi_table = pd.read_table(self.voi, sep=',')
        elif os.path.splitext(self.voi)[1] == '.tsv':
            voi_table = pd.read_table(self.voi, sep='\t')
        #Concatenate gene, ref, alt and position to create unique index string
        voi_table['Variant'] = voi_table['Gene'] + ':' + voi_table['RefAA'] + \
            voi_table['AAPos'].astype(str) + voi_table['AltAA']
        #Concatenate ref, alt and position to create unique SNP identifier
        voi_table['SNP'] = voi_table['RefAA']+voi_table['AAPos'].astype(str) + \
            voi_table['AltAA']
        voi_table.set_index(['Variant'], inplace=True)
        return(voi_table)

    def getGeneStats(self, bam_path):
        """Given a BAM file and BED file, determine the per base coverage for the
        each chromosome in the BED file."""
        bam_file = pysam.AlignmentFile(bam_path, 'rb')
        bed = Bed(self.bed)
        gene_stats = OrderedDict()
        for bed_rec in bed.getExonTable():
            try:
                gene_stats[bed_rec.chrom] += list(np.sum(
                    bam_file.count_coverage(bed_rec.chrom,
                    bed_rec.start, bed_rec.stop+1), axis=0))
            except KeyError:
                gene_stats[bed_rec.chrom] = list(np.sum(
                    bam_file.count_coverage(bed_rec.chrom, bed_rec.start,
                    bed_rec.stop+1), axis=0))
        return(gene_stats)

    def getExpCoverage(self):
        """Given a study, return a dictionary of dictionary storing the per base
        coverage per gene per sample."""
        bam_files = glob.glob('{0}/*/alignments/output_FM_SR_DD_RG.bam'.format(
            self.out_path))
        barcode = re.compile('_[ATGC]*-[ATGC]*')
        sample_gstat = OrderedDict()
        for files in bam_files:
            sample_dir = os.path.basename(os.path.dirname(files))
            sample = barcode.split(sample_dir)[0]
            sample_gstat[sample] = self.getGeneStats(files)
        return(sample_gstat)

    def checkDepthPass(self):
        """Given a study, return a dictionary of dictionary describing whether
        the gene passed depth threshold per gene per study."""
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
        """Given a study output path, generate a pandas dataframe containing
        all intronic variants. The module gets a list of variants from the
        output folder. A dictionary of list is initialzed with the fields that
        will be recorded in the intron summary tables. The features of the
        intron summary tables are listed below:
            1. Chrom: Chromosome where the variant is found.
            2. Gene: Gene boundary within which the variant is found.
            3. Pos: Variant position.
            4. Qual: Phred based quality score of variant call.
            5. Ref: Reference base call for the variant call.
            6. Alt: Alternate(variant) base call for the variant call.
            7. Exon: Exon boundary within which the variant is found. This
                     value is set to "Intron" for the intron table.
            8. AAPos: Amino acid position for the variant call. This value is
                      set to np.nan value for the intron table.
            9. RefCodon: Reference triplet codon for the variant call. This
                         value is set to "NA" for the intron table.
           10. AltCodon: Alternate(variant) base call for the variant call. This
                         value is set to "NA" for the intron table.
           11. RefAA: Reference amino acid for the variant call. This value is
                      set to "NA" for the intron table.
           12. AltAA: Alternate(variant) base call for the variant call. This
                      value is set to "NA" for the intron table.
           13. DP: Depth for the codon to which the variant call belong.
           14. AF: Allele frequency for the variant call.
           15. Conf: Number of variant callers that identified the variant call.
        The table is indexed by the following fields:
            1. Sample: Sample name
            2. Variant: A string containing the variant information in
                        <Chrom>:<Ref><Pos><Alt> for the intron table. """

        vcf_files = glob.glob('{0}/*/*_variants_merged_annotated.vcf'.format(
            self.out_path))
        vcf_dict = {'Gene': [], 'Pos': [], 'Qual': [], 'Ref': [], 'Alt': [],
            'AAPos': [], 'AltCodon': [], 'RefCodon': [], 'RefAA': [],
            'AltAA': [], 'DP': [], 'AF': [], 'Confidence': [], 'Exon': [],
            'Chrom': [], 'Sources': []}
        vcf_var = list()
        vcf_sample = list()
        for files in vcf_files:
            vcf = Reader(files, self.fasta)
            vcf.readheader()
            vcf_file = vcf.readvcf()
            barcode = re.compile('_[ATGC]*-[ATGC]*')
            for var in vcf_file:
                #If variant call is annotated as an intronic call
                #push it into the intronic variant dictionary
                sample = barcode.split(list(var.Samples.keys())[0])[0]
                variant = '{0}:{1}{2}{3}'.format(var.CHROM, var.REF[0],
                                                var.POS, var.ALT[0])
                if 'Exon' not in var.INFO:
                    continue
                elif var.INFO['Exon'][0] == 'Intron':
                    vcf_dict['Chrom'].append(var.CHROM)
                    vcf_dict['Gene'].append(var.CHROM)
                    vcf_dict['Pos'].append(var.POS)
                    vcf_dict['Qual'].append(var.QUAL)
                    vcf_dict['Ref'].append(var.REF[0])
                    vcf_dict['Alt'].append(var.ALT[0])
                    vcf_dict['Exon'].append('Intron')
                    vcf_dict['AAPos'].append(np.nan)
                    vcf_dict['RefCodon'].append('NA')
                    vcf_dict['AltCodon'].append('NA')
                    vcf_dict['RefAA'].append('NA')
                    vcf_dict['AltAA'].append('NA')
                    try:
                        vcf_dict['DP'].append(var.INFO['DP'][0])
                    except KeyError:
                        vcf_dict['DP'].append(0)
                    vcf_dict['AF'].append(float(var.INFO['Freq'][0])*100)
                    vcf_dict['Confidence'].append(int(var.INFO['Confidence'][0]))
                    vcf_dict['Sources'].append(','.join(var.INFO['Sources']))
                    vcf_var.append(variant)
                    vcf_sample.append(sample)
        vcf_index = [np.array(vcf_sample), np.array(vcf_var)]
        vcf_df = pd.DataFrame(vcf_dict, index=vcf_index)
        vcf_df.index.names = ['Sample', 'Variant']
        return(vcf_df)

    def getVarTables(self):
        """Given a study output path, generate a pandas dataframe containing
        all known and novel exonic variants. The module gets a list of VCF files
        from the output folder. A dictionary of list is initialzed with the
        fields that will be recorded in the intron summary tables. If a variant
        from the list of variants of intrest is not found in the sample, an
        empty record is added to the dataframe with np.nan or "NA" values for
        all field expect the variant descriptor fields. The features
        of the exonic summary tables are listed below:
            1. Chrom: Chromosome where the variant is found.
            2. Gene: Gene boundary within which the variant is found.
            3. Pos: Variant position.
            4. Qual: Phred based quality score of variant call.
            5. Ref: Reference base call for the variant call.
            6. Alt: Alternate(variant) base call for the variant call.
            7. Exon: Exon boundary within which the variant is found. This
                     value is set to exon number for the known and novel exonic
                     table.
            8. AAPos: Amino acid position for the variant call.
            9. RefCodon: Reference triplet codon for the variant call.
           10. AltCodon: Alternate(variant) base call for the variant call.
           11. RefAA: Reference amino acid for the variant call.
           12. AltAA: Alternate(variant) base call for the variant call.
           13. DP: Depth for the codon to which the variant call belong.
           14. AF: Allele frequency for the variant call.
           15. Conf: Number of variant callers that identified the variant call.
        The table is indexed by the following fields:
            1. Sample: Sample name
            2. Variant: A string containing the variant information in
                        <Chrom>:<RefAA><AAPos><AltAA> for the intron table. """

        voi_table = self.getVarOfInt()
        vcf_files = glob.glob('{0}/*/*_variants_merged_annotated.vcf'.format(
                                                                 self.out_path))
        vcf_df = pd.DataFrame()
        vcf_dict = {'Gene' : [], 'Pos' : [], 'Qual' : [], 'Ref' : [],
            'Alt' : [], 'AAPos' : [], 'RefCodon' : [], 'AltCodon' : [],
            'RefAA' : [], 'AltAA' : [], 'DP' : [], 'AF' : [], 'Confidence': [],
            'Exon' : [], 'Chrom' : [], 'Sources' : []}
        vcf_var = list()
        vcf_sample = list()
        vcf_gene = list()
        var_sample = list()
        voi_df = self.getVarOfInt()
        for files in vcf_files:
            vcf = Reader(files, self.fasta)
            vcf.readheader()
            vcf_file = vcf.readvcf()
            barcode = re.compile('_[ATGC]*-[ATGC]*')
            count = 0
            for var in vcf_file:
                sample = barcode.split(list(var.Samples.keys())[0])[0]
                if 'Gene' not in var.INFO:
                    continue
                variant = '{0}:{1}{2}{3}'.format(var.INFO['Gene'][0],
                                                var.INFO['RefAA'][0],
                                                var.INFO['AAPos'][0],
                                                var.INFO['AltAA'][0])
                sample_variant = '{0}{1}'.format(sample, variant)
                if (var.INFO['Gene'][0] == None or
                        var.INFO['RefAA'][0] == None or
                        var.INFO['AAPos'][0] == None or
                        var.INFO['AltAA'][0] == None):
                    continue
                if variant in voi_df.index:
                    count += 1
                if sample_variant in var_sample:
                    index = 0
                    for pos in range(len(vcf_var)):
                        if sample == vcf_sample[pos] and variant == vcf_var[pos]:
                            index = pos
                    vcf_dict['Ref'][index] = '{0},{1}'.format(
                                             vcf_dict['Ref'][index], var.REF[0])
                    vcf_dict['Alt'][index] = '{0},{1}'.format(
                                             vcf_dict['Alt'][index], var.ALT[0])
                else:
                    vcf_dict['Chrom'].append(var.CHROM)
                    vcf_dict['Gene'].append(var.INFO['Gene'][0])
                    vcf_dict['Pos'].append(var.POS)
                    vcf_dict['Qual'].append(var.QUAL)
                    vcf_dict['Ref'].append(var.REF[0])
                    vcf_dict['Alt'].append(var.ALT[0])
                    vcf_dict['Exon'].append(var.INFO['Exon'][0])
                    vcf_dict['AAPos'].append(int(var.INFO['AAPos'][0]))
                    vcf_dict['RefCodon'].append(var.INFO['RefCodon'][0])
                    vcf_dict['AltCodon'].append(var.INFO['AltCodon'][0])
                    vcf_dict['RefAA'].append(var.INFO['RefAA'][0])
                    vcf_dict['AltAA'].append(var.INFO['AltAA'][0])
                    vcf_dict['DP'].append(var.INFO['DP'][0])
                    vcf_dict['AF'].append(float(var.INFO['Freq'][0]) * 100)
                    vcf_dict['Confidence'].append(int(var.INFO['Confidence'][0]))
                    vcf_dict['Sources'].append(','.join(var.INFO['Sources']))
                    vcf_gene.append(var.CHROM)
                    vcf_var.append(variant)
                    vcf_sample.append(sample)
                    try:
                        if var.INFO['DP'][0] > 0:
                            var_sample.append(sample_variant)
                        else:
                            var_sample.append('{0}{1}:{2}{3}NA'.format(sample,
                                                            var.INFO['Gene'][0],
                                                           var.INFO['RefAA'][0],
                                                           var.INFO['AAPos'][0],
                                                           var.INFO['AltAA'][0]))
                    except TypeError:
                        if var.INFO['DP'][0] is None:
                           var_sample.append('{0}{1}:{2}{3}NA'.format(sample,
                                                             var.INFO['Gene'][0],
                                                            var.INFO['RefAA'][0],
                                                            var.INFO['AAPos'][0],
                                                            var.INFO['AltAA'][0]))
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
                    vcf_dict['Confidence'].append(3)
                    vcf_dict['Sources'].append('GATK,Freebayes,Samtools')
                    vcf_var.append(variants)
                    vcf_sample.append(sample)
        vcf_index = [np.array(vcf_sample), np.array(vcf_var)]
        vcf_df = pd.DataFrame(vcf_dict, index=vcf_index)
        vcf_df.index.names = ['Sample', 'Variant']
        return(vcf_df)

    def getRepSnps(self):
        """Returns a table containing sample wise breakdown about the presence
        or absence of variants of interest in the study. If there are samples
        with no calls for any particular variants of interest, if the sample has
        coverage at the variant location and no variant call, record is labeled
        as a "WT" call in the FinalCall field. If there is a variant call at any
        of the locations, the record is labeled as "SNP" call in the FinalCall
        field. This module internally calls the following modules:
            1. getVarTables()
            2. getVarOfInt()

        The module produces a pandas dataframe with the following fields:
          1. Sample : Sample name
          2. Variant : Variant of interest
          3. Chrom : Chromosome
          4. Gene : Gene name
          5. SNP: SNP of interest
          6. FinalCall : Whether the sample had the SNP or it was a wild type call (WT)
          7. Ref : Reference genomic base
          8. Alt : Alternate(Variant) genomic base
          9. Pos : Genomic position of variant
          10. Qual : Phred based quality score
          11. RefCodon : Reference codon at variant of interest
          12. RefAA : Reference amino acid at variant of interest
          13. AltCodon : Alternate codon at variant of interest
          14. AltAA : Alternate amino acid at variant of interest
          15. AAPos : Amino acid position of interest
          16. Exon : Exon number where the variant is found
          17. AF : Allele frequency of variant of interest
          18. DP : Depth of coverage of variant of interest
          19. Conf : Number of variant callers that called the variant of interest
        """
        #Get table of exonic variants and variants of interest
        pd.set_option('display.max_columns', None)
        exp_df = self.getVarTables()
        voi_df = self.getVarOfInt()
        exp_voi = pd.DataFrame()
        if voi_df is None:
            return(None)
        for sample, var_df in exp_df.groupby(level=0):
            sam_index = list()
            var_df = var_df.reset_index(level=0)
            var_voi = var_df.merge(voi_df, how='right', left_index=True,
                right_index=True)
            #Create a list of length equal to number of variants of interest
            #containing the sample name, to ensure that the final table has
            #exactly the same number of records per sample
            sam_index = [sample] * len(var_voi)
            var_index = [np.array(sam_index), np.array(var_voi.index)]
            var_voi.set_index(var_index, inplace=True)
            var_voi.index.names = ['Sample', 'Variant']
            exp_voi = exp_voi.append(var_voi)
            
        exp_voi['FinalCall'] = exp_voi['SNP']
        #Regex to check if the variant description field is in the correct format
        var_regex = (r'(?P<RefAA>[DTSEPGACVMILYFHKRWQN])'
                     r'(?P<AAPos>\d+)(?P<AltAA>[DTSEPGACVMILYFHKRWQN])')
        for index, series in exp_voi.iterrows():
            if pd.isnull(series['DP']) or series['DP'] == 0:
                exp_voi.at[index, 'FinalCall'] = 'WT'
                exp_voi.at[index, 'Confidence'] = 3
                exp_voi.at[index, 'Sources'] = 'GATK,Samtools,Freebayes'
            elif pd.isnull(series['Alt']):
                var_reg = re.match(var_regex, series['SNP'])
                exp_voi.at[index, 'FinalCall'] = '{0}{1}{0}'.format(
                                                         var_reg.group('RefAA'),
                                                         var_reg.group('AAPos'))
            exp_voi.at[index, 'AAPos'] = series['AAPos_x'] and series['AAPos_y']
            exp_voi.at[index, 'AltAA'] = series['AltAA_x'] and series['AltAA_y']
            exp_voi.at[index, 'Chrom'] = series['Chrom_x'] and series['Chrom_y']
            exp_voi.at[index, 'Gene'] = series['Gene_x'] and series['Gene_y']
            exp_voi.at[index, 'RefAA'] = series['RefAA_x'] and series['RefAA_y']
        #try:
        exp_voi.drop(['AAPos_x', 'AAPos_y', 'AltAA_x', 'AltAA_y',
            'Chrom_x', 'Chrom_y', 'Gene_x', 'Gene_y', 'RefAA_x', 'RefAA_y'],
            inplace=True, axis=1)
        exp_voi = exp_voi[['Chrom', 'Gene', 'SNP', 'FinalCall', 'Ref', 'Alt',
            'Pos', 'Qual', 'RefCodon', 'RefAA', 'AltCodon', 'AltAA', 'AAPos',
            'Exon', 'AF', 'DP', 'Confidence', 'Sources']]
        #Sorting to remove performance warning
        exp_voi.sort_index(inplace=True)
        return(exp_voi)

    def getNovSnps(self):
        """Returns a table containing sample wise breakdown about the presence
        of novel variants in the study. This table will contain a FinalCall
        field, since only variants are reported and not wild type calls. This
        module internally calls the following modules:
            1. getVarTables()
            2. getVarOfInt()
        The module produces a pandas dataframe with the following fields:
          1. Sample : Sample name
          2. Variant : Variant of interest
          3. Chrom : Chromosome
          4. Gene : Gene name
          5. Ref : Reference genomic base
          6. Alt : Alternate(Variant) genomic base
          7. Pos : Genomic position of variant
          8. Qual : Phred based quality score
          9. RefCodon : Reference codon at variant of interest
          10. RefAA : Reference amino acid at variant of interest
          11. AltCodon : Alternate codon at variant of interest
          12. AltAA : Alternate amino acid at variant of interest
          13. AAPos : Amino acid position of interest
          14. Exon : Exon number where the variant is found
          15. AF : Allele frequency of variant of interest
          16. DP : Depth of coverage of variant of interest
          17. Conf : Number of variant callers that called the variant of interest
        """

        exp_df = self.getVarTables()
        voi_df = self.getVarOfInt()
        exp_nov = pd.DataFrame()
        for sample, var_df  in exp_df.groupby(level=0):
            sam_index = list()
            var_df = var_df.reset_index(level=0)
            if voi_df is None:
                var_nov = var_df
            else:
                var_nov = var_df[~var_df.index.isin(voi_df.index)]
            #Create a list of length equal to number of novel variants for that
            #sample, containing the sample name. This enables the table to be
            #indexed by sample name.
            sam_index = [sample] * len(var_nov)
            var_index = [np.array(sam_index), np.array(var_nov.index)]
            var_nov.set_index(var_index, inplace=True)
            var_nov.index.names = ['Sample', 'Variant']
            exp_nov = exp_nov.append(var_nov)
        exp_nov = exp_nov[exp_nov.Confidence >= 2]
        exp_nov = exp_nov[['Chrom', 'Gene', 'Ref', 'Alt', 'Pos', 'Qual',
            'RefCodon', 'RefAA', 'AltCodon', 'AltAA', 'AAPos', 'Exon', 'AF',
            'DP', 'Confidence', 'Sources']]

        return(exp_nov)

    def getBamStat(self, bamfile, chrom, start, stop):
        """Given BAM file path, chromosome name, start and stop, return
        average coverage for the chromosome."""
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
        """Given a data frame containing known variant calls, edit the depth to
        reflect average codon for variant of interest."""
        depth_list = list()
        for row, value in var_df.iterrows():
            bamfile = glob.glob('{0}/{1}*/alignments/output_FM_SR_DD_RG.bam'.format(
                                                      self.out_path, row[0]))[0]
            nuc_pos = self.getBaseRange(value.Chrom, value.Gene,
                        value.AAPos)
            if nuc_pos == None:
                nuc_pos = range(int(value.Pos) -1, int(value.Pos) + 2)
            depth = self.getBamStat(bamfile, value.Chrom, nuc_pos.start,
                nuc_pos.stop)
            depth_list.append(depth)
                   #np.log10(depth+1))
        var_df['DP'] = pd.Series(depth_list, index=var_df.index)
        return(var_df)

    def toJSON(self):
        """Create a summary compressed JSON file for all variants in all samples
        in the study. The information recorded in the JSON format include:
            1. Study name
            2. Number of samples in study
            3. Date of analysis
            4. Sample name:
                i. Variant description string in <Gene><RefAA><AAPos><AltAA>
                   format:
                   a. Reference amino acid
                   b. Alternate amino acid
                   c. Variant call position
                   d. Final call field, indicating the wild type or variant call
                      based of BAM file coverage
                   e. Allele frequency for the variant call
                   f. Depth for the variant call
                   g. Confidence for the vairant call based on the number of
                      variant callers that identified the variant
        """

        known_snps = self.getRepSnps()
        novel_snps = self.getNovSnps()
        known_snps = self.getDepthStats(known_snps)
        novel_snps = self.getDepthStats(novel_snps)
        study = os.path.basename(self.out_path)
        sampleCount = len(glob.glob(
                 '{0}/*/*_variants_merged_annotated.vcf'.format(self.out_path)))
        analysisTime = datetime.now()
        analysisDate = analysisTime.strftime('%m-%d-%y')
        jsonFile = open('{0}/Study_variants.json'.format(self.out_path), 'w')
        json_dict = OrderedDict({'Study': study,
            'Number of samples': sampleCount, 'Date': analysisDate})

        if known_snps is None:
            self.logger.debug('No vairants of interest provided; generating only novel SNP summary')
        else:
            for index, record in known_snps.iterrows():
                sample_info = self.config[index[0]]
                if 'Sample' not in json_dict:
                    json_dict['Sample'] = {}
                if index[0] not in json_dict['Sample']:
                    json_dict['Sample'][index[0]] = {}
                if 'VariantCalls' not in json_dict['Sample'][index[0]]:
                    json_dict['Sample'][index[0]]['Year'] = sample_info.year
                    json_dict['Sample'][index[0]]['Country'] = sample_info.country
                    json_dict['Sample'][index[0]]['Site'] = sample_info.site
                    json_dict['Sample'][index[0]]['TreatmentDay'] = sample_info.treatmentDay
                    json_dict['Sample'][index[0]]['SampleID'] = sample_info.iD
                    json_dict['Sample'][index[0]]['Genus&Species'] = sample_info.genus
                    json_dict['Sample'][index[0]]['SampleType'] = sample_info.type
                    json_dict['Sample'][index[0]]['Markers'] = ','.join(sample_info.markers.tolist())
                    json_dict['Sample'][index[0]]['Replicate'] = sample_info.replicate
                    json_dict['Sample'][index[0]]['VariantCalls'] = {}
                json_dict['Sample'][index[0]]['VariantCalls'][index[1]] = {
                    'Ref' : record.RefAA, 'Pos': record.AAPos,
                    'Alt': record.AltAA , 'Call' :  record.FinalCall,
                    'AF' : record.AF, 'DP': record.DP, 'Confidence': record.Confidence,
                    'Status' : 'Known', 'Sources': record.Sources}

        for index, record in novel_snps.iterrows():
            sample_info = self.config[index[0]]
            if 'Sample' not in json_dict:
                json_dict['Sample'] = OrderedDict()
            if index[0] not in json_dict['Sample']:
                json_dict['Sample'][index[0]] = OrderedDict()
            if 'VariantCalls' not in json_dict['Sample'][index[0]]:
                json_dict['Sample'][index[0]]['Year'] = sample_info.year
                json_dict['Sample'][index[0]]['Country'] = sample_info.country
                json_dict['Sample'][index[0]]['Site'] = sample_info.site
                json_dict['Sample'][index[0]]['TreatmentDay'] = sample_info.treatmentDay
                json_dict['Sample'][index[0]]['SampleID'] = sample_info.iD
                json_dict['Sample'][index[0]]['Genus&Species'] = sample_info.genus
                json_dict['Sample'][index[0]]['SampleType'] = sample_info.type
                json_dict['Sample'][index[0]]['Markers'] = ','.join(sample_info.markers.tolist())
                json_dict['Sample'][index[0]]['Replicate'] = sample_info.replicate
                json_dict['Sample'][index[0]]['VariantCalls'] = OrderedDict()
            calls = '{0}{1}{2}'.format(record.RefAA, record.AAPos, record.AltAA)
            json_dict['Sample'][index[0]]['VariantCalls'][index[1]] = {
                'Ref' : record.RefAA, 'Pos': record.AAPos,
                'Alt': record.AltAA , 'Call' :  'NA',
                'AF' : record.AF, 'DP': record.DP, 'Confidence': record.Confidence,
                'Status' : 'Novel', 'Sources': record.Sources}

        json.dump(json_dict, jsonFile, indent=4)
        jsonFile.close()
        return('{0}/Study_variants.json'.format(self.out_path))

    def toCSV(self, var_type):
        """Create CSV files from the DataFrames and generate the allele
        frequency and depth views for known and novel files"""
        #Sumarize variants of intrest
        rep_dir = '{0}/Reports'.format(self.out_path)
        if not os.path.exists(rep_dir):
            os.mkdir(rep_dir)
        if var_type == 'known':
            var_df = self.getRepSnps()
            if var_df is None:
                self.logger.debug('No known variants provided')
            else:
                var_df = self.getDepthStats(var_df)
                var_df = var_df.reset_index(level=1)
                var_df['Test'] = var_df.index
                var_df.drop_duplicates(['Test', 'Variant'], inplace=True)
                var_df.drop('Test', axis=1, inplace=True)
                out_file = '{0}/Study_known_variants.csv'.format(rep_dir)
        elif var_type == 'novel':
            var_df = self.getNovSnps()
            var_df = self.getDepthStats(var_df)
            var_df = var_df.reset_index(level=1)
            out_file = '{0}/Study_novel_exonic_variants.csv'.format(
                rep_dir)
            sout_file = '{0}/Study_novel_intronic_variants.csv'.format(
                rep_dir)
            exp_intron = self.getIntronTables()
            exp_intron = exp_intron.reset_index()
            intron_key = ['Gene_name', 'RefAA_sym', 'AAPos_sort', 'AltAA_sym']
            intron_regex = ('(?P<Gene_name>[a-zA-Z0-9]+):'
                            '(?P<RefAA_sym>[a-zA-Z]?)(?P<AAPos_sort>[0-9]+)'
                            '(?P<AltAA_sym>[a-zA-Z]?)')
            try:
                exp_intron[intron_key] = exp_intron['Variant'].str.extract(
                    intron_regex, expand=True)
                exp_intron['AAPos_sort'] = pd.to_numeric(exp_intron['AAPos_sort'])
                exp_intron.sort_values(['Sample', 'Gene_name', 'AAPos_sort'],
                    inplace=True)
                exp_intron.drop(labels=['Gene_name', 'RefAA_sym', 'AAPos_sort',
                    'AltAA_sym' ], axis=1, inplace=True)
                exp_intron.sort_index().reset_index(drop=True).to_csv(sout_file,
                     index=False)
            except AttributeError:
                self.logger.debug('No novel intronic variants found')

        if var_df is not None:
            var_key = ['Gene_name', 'RefAA_sym', 'AAPos_sort', 'AltAA_sym']
            var_regex = ('(?P<Gene_name>[a-zA-Z0-9]+):'
                         '(?P<RefAA_sym>[a-zA-Z]?)(?P<AAPos_sort>[0-9]+)'
                         '(?P<AltAA_sym>[a-zA-Z]?)')
            try:
                var_df[var_key] = var_df['Variant'].str.extract(var_regex, expand=True)
                var_df['Sample_name'] = var_df.index
                var_df['AAPos_sort'] = pd.to_numeric(var_df['AAPos_sort'])
                var_df.sort_values(['Sample_name', 'Gene_name', 'AAPos_sort'],
                    inplace=True)
                var_df.drop(labels=['Sample_name', 'Gene_name', 'RefAA_sym',
                    'AAPos_sort', 'AltAA_sym'], axis=1, inplace=True)
                var_df.to_csv(out_file)
                exp_af = var_df.pivot(var_df.index, 'Variant')['AF'].transpose()
                exp_af['Variant'] = exp_af.index
                exp_af[var_key] = exp_af['Variant'].str.extract(var_regex, expand=True)
                exp_af['AAPos_sort'] = pd.to_numeric(exp_af['AAPos_sort'])
                exp_af.sort_values(['Gene_name', 'AAPos_sort'], inplace=True)
                exp_af.drop(labels=['Variant', 'Gene_name', 'RefAA_sym', 'AAPos_sort',
                          'AltAA_sym'], axis=1, inplace=True)
                af_mask = exp_af.isnull()
                exp_af.to_csv('{0}/Study_{1}_variants_allele_frequency.csv'.format(
                                                           rep_dir, var_type))

                exp_dp = var_df.pivot(var_df.index, 'Variant')['DP'].transpose()
                exp_dp['Variant'] = exp_dp.index
                exp_dp[var_key] = exp_dp['Variant'].str.extract(var_regex, expand=True)
                exp_dp['AAPos_sort'] = pd.to_numeric(exp_dp['AAPos_sort'])
                exp_dp.sort_values(['Gene_name', 'AAPos_sort'], inplace=True)
                exp_dp.drop(labels=['Variant', 'Gene_name', 'RefAA_sym', 'AAPos_sort',
                          'AltAA_sym'], axis=1, inplace=True)
                dp_mask = exp_dp.isnull()
                exp_dp.to_csv('{0}/Study_{1}_variants_depth.csv'.format(rep_dir,
                                                                          var_type))
            except AttributeError:
                self.logger.debug('No novel exonic variants found') 

    def getSummary(self):
        """Generate CSV tables, JSON and figures for a given study"""
        #Write Known and novel variants to files
        self.toCSV('known')
        self.toCSV('novel')
        fig_path = '{0}/Figures'.format(self.out_path)
        if not os.path.exists(fig_path):
            os.mkdir(fig_path)
        # Plot using Rscript
        self.logger.info('Plotting Depth Per SNP')
        dcmd = ['Rscript',
            '{0}/Rscripts/DepthPerReportSNP.R'.format(self.summary_path), '-i',
            '{0}/Reports/Study_known_variants_depth.csv'.format(self.out_path),
            '-o', '{0}/Study_depth.pdf'.format(fig_path)]
        drun = subprocess.Popen(dcmd, shell=False,
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        drun.wait()
        if drun.returncode != 0:
            self.logger.error('Failed to execute DepthPerReportSNP.R')
            self.logger.error(' '.join(dcmd))

        self.logger.info('Plotting Reportable SNPs Frequency')
        acmd = ['Rscript',
            '{0}/Rscripts/reportableSNPsFreq.R'.format(self.summary_path), '-i',
            '{0}/Reports/Study_known_variants.csv'.format(self.out_path), '-r',
            '{0}'.format(self.voi),
            '-o', '{0}/'.format(fig_path)]
        arun = subprocess.Popen(acmd, shell=False,
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        arun.wait()
        if arun.returncode != 0:
            self.logger.error('Failed to execute reportableSNPsFreq.R')
            self.logger.error(' '.join(acmd))
        self.logger.info('Plotting Novel Exonic Non-Synonymous SNPs')
        nenscmd = ['Rscript',
            '{0}/Rscripts/NovelExonicNonSynSNPs.R'.format(self.summary_path),
            '-i', '{0}/Reports/Study_novel_exonic_variants.csv'.format(self.out_path),
            '-o', '{0}/'.format(fig_path)]
        nensrun = subprocess.Popen(nenscmd, shell=False,
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        nensrun.wait()
        if nensrun.returncode != 0:
            self.logger.warning('Failed to execute NovelExonicNonSynSNPs.R; This could be because no novel non-synonymous variants were found, check Study_novel_exonic_variants.csv')
            self.logger.warning(' '.join(nenscmd))

        self.logger.info('Plotting Novel Exonic Synonymous SNPs')
        nescmd = ['Rscript',
            '{0}/Rscripts/NovelExonicSynSNPs.R'.format(self.summary_path), '-i',
            '{0}/Reports/Study_novel_exonic_variants.csv'.format(self.out_path),
            '-o', '{0}/'.format(fig_path)]
        nesrun = subprocess.Popen(nescmd, shell=False,
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        nesrun.wait()
        if nesrun.returncode != 0:
            self.logger.warning('Failed to execute NovelExonicSynSNPs.R; This could be because no novel synonymous variants were found, check Study_novel_exonic_variants.csv')
            self.logger.warning(' '.join(acmd))

        self.logger.info('Plotting Novel Intronic SNPs')
        nicmd = ['Rscript',
            '{0}/Rscripts/NovelIntronicSNPs.R'.format(self.summary_path), '-i',
            '{0}/Reports/Study_novel_intronic_variants.csv'.format(self.out_path),
                '-o', '{0}/'.format(fig_path)]
        nirun = subprocess.Popen(nicmd, shell=False,
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        nirun.wait()
        if nirun.returncode != 0:
            self.logger.warning('Failed to execute NovelIntronicSNPs.R; This could be because no novel intronic variants were found, check Study_novel_intronic_variants.csv')
            self.logger.warning(' '.join(nicmd))

        #os.remove('{0}/Reportable_SNPs_Report.csv'.format(out_dir))
        if os.path.exists('{0}/novel_SNPs_exonic_syn.csv'.format(fig_path)):
            os.remove('{0}/novel_SNPs_exonic_syn.csv'.format(fig_path))
        if os.path.exists('{0}/novel_SNPs_intronic.csv'.format(fig_path)):
            os.remove('{0}/novel_SNPs_intronic.csv'.format(fig_path))
        if os.path.exists('{0}/novel_SNPs_exonic_nonsyn.csv'.format(fig_path)):
            os.remove('{0}/novel_SNPs_exonic_nonsyn.csv'.format(fig_path))
        if os.path.exists('{0}/Study_novel_exonic_variants_filtered.csv'.format(fig_path)):
            os.remove('{0}/Study_novel_exonic_variants_filtered.csv'.format(
                      fig_path))
        if os.path.exists('{0}/Study_novel_intronic_variants_filtered.csv'.format(fig_path)):
            os.remove('{0}/Study_novel_intronic_variants_filtered.csv'.format(
                      fig_path))

    def getVarStats(self, vcf_file):
        vcf_file = Reader(vcf_file, self.fasta)
        vcf_file.readheader()
        vcf_reader = vcf_file.readvcf()
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
            try:
                total += 1
                if variant.INFO['Confidence'][0] >= 2:
                    verfied += 1
                if variant.INFO['Exon'][0] == 'Intron':
                    intronic += 1
                else:
                    exonic += 1
                    if variant.INFO['RefAA'][0] == variant.INFO['AltAA'][0]:
                        syn += 1
                    else:
                        nsyn += 1
                    if '{0}{1}'.format(variant.REF[0], str(variant.ALT[0])) in trasition:
                        trans += 1
                    else:
                        tranv += 1
            except KeyError:
                continue
        return(total, verfied, exonic, intronic, syn, nsyn, trans, tranv)
