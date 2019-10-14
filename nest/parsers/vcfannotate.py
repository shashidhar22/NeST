import os
import math
import pysam 
import logging
import pandas as pd
import numpy as np
from collections import namedtuple
from collections import OrderedDict
from nest.parsers.bed import Bed
from nest.parsers.fasta import Fasta
from nest.parsers.vcfReader import Reader 
from nest.parsers.vcfwriter import Writer

class Annotate:
    def __init__(self):
        self.logger = logging.getLogger('NeST.VcfAnnotate')
        #self.vcf = vcf_object
        return

    def addInfoHeader(self, name, description=None, version=None,
                source=None, type='String', number=1):
        info_field = namedtuple('Info', ['id', 'number', 'type',
                                            'description',
                                            'version', 'source'])
        info_field.id = name
        info_field.number = number
        info_field.type = type
        info_field.description = description
        info_field.version = version
        info_field.source = source
        return(info_field)

    def getAA(self, codon):
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
            'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
            }
        amino_acid = codontable[codon]
        return(amino_acid)

    def peek(self, vcf_reader):
        curr_pos = vcf_reader.tell()
        curr_line = vcf_reader.readline()
        vcf_reader.seek(curr_pos)
        return(curr_line)

    def getTriplet(self, bed, fasta_path, uid, cds_pos):
        exon_table = bed.getExonTable()
        exon = ''
        strand = '+'
        for records in exon_table:
            if uid in range(records.uidStart, records.uidStop):
                exon = records.exon
                strand = records.strand
        coding_fasta = bed.getCodingFasta(fasta_path)
        for seq in coding_fasta:
            if uid in range(seq.fid, seq.fid + seq.length):
                #print(uid, cds_pos, seq.fid, seq.fid+seq.length+1, seq.chrom)
                codon_pos = cds_pos % 3
                if codon_pos == 0:
                    aa_pos = int(cds_pos/3)
                    #Account for strand and return corresponding co
                    if strand == '+':
                        return(seq.seq[cds_pos-3: cds_pos], aa_pos, codon_pos,
                                exon, seq.chrom, strand)
                    else:
                        return(seq.seq[::-1][cds_pos-3: cds_pos], aa_pos, codon_pos,
                                exon, seq.chrom, strand)
                elif codon_pos == 1:
                    aa_pos = int(cds_pos/3) + 1
                    if strand == '+':
                        return(seq.seq[cds_pos-1: cds_pos+2], aa_pos, codon_pos,
                                exon, seq.chrom, strand)
                    else:
                        return(seq.seq[::-1][cds_pos-1: cds_pos+2], aa_pos, codon_pos,
                                exon, seq.chrom, strand)
                elif codon_pos == 2:
                    aa_pos = int(cds_pos/3) + 1
                    if strand == '+':
                        return(seq.seq[cds_pos-2: cds_pos+1], aa_pos, codon_pos,
                                exon, seq.chrom, strand)
                    else:
                        return(seq.seq[::-1][cds_pos-2: cds_pos+1], aa_pos, codon_pos,
                                exon, seq.chrom, strand)

    def getModFasta(self, fasta_path, vcf_path):
        vcf = Reader(vcf_path, fasta_path)
        vcf.readheader()
        vcf_reader = vcf.readvcf()
        fasta = Fasta(fasta_path)
        fasta_reader = fasta.read()
        for fasta_rec in fasta_reader:
            vcf = Reader(vcf_path, fasta_path)
            vcf.readheader()
            vcf_reader = vcf.readvcf()
            header = fasta_rec.header
            seq = fasta_rec.seq
            fid = fasta_rec.fid
            length = fasta_rec.length
            for vcf_rec in vcf_reader:
                if vcf_rec.UID in range(fasta_rec.fid,
                                fasta_rec.fid + fasta_rec.length ):
                    if len(vcf_rec.REF[0]) > 1 or len(vcf_rec.ALT[0]) > 1:
                        continue
                    else:
                        seq = seq[:vcf_rec.POS-1] + vcf_rec.ALT[0] +\
                                                    seq[vcf_rec.POS:]
            mod_seq = namedtuple('fastaRec', ['header', 'seq', 'fid',
                                            'length'])
            record = mod_seq(header, seq, fid, length)
            yield record

    def getAltTriplet(self, bed, mod_fasta, vcf_rec, codon_pos, strand):
        fasta = Fasta(mod_fasta).read()
        comp = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
        if strand == '+':
            codon_pos = codon_pos
        else:
            codon_pos = 3 - codon_pos
            if codon_pos == 3:
                codon_pos = 0
        for rec in fasta:
            if rec.header == vcf_rec.CHROM:
                if codon_pos == 0:
                    if strand == '+':
                        return(rec.seq[vcf_rec.POS-3 : vcf_rec.POS])
                    else:
                        codon = rec.seq[vcf_rec.POS-3 : vcf_rec.POS][::-1]
                        codon = ''.join([comp[nuc] for nuc in codon])
                        return(codon)
                elif codon_pos == 1:
                    if strand == '+':
                        return(rec.seq[vcf_rec.POS-1 : vcf_rec.POS+2])
                    else:
                        codon = rec.seq[vcf_rec.POS-1 : vcf_rec.POS+2][::-1]
                        codon = ''.join([comp[nuc] for nuc in codon])
                        return(codon)
                elif codon_pos == 2:
                    if strand == '+':
                        return(rec.seq[vcf_rec.POS-2 : vcf_rec.POS+1])
                    else:
                        codon = rec.seq[vcf_rec.POS-2 : vcf_rec.POS+1][::-1]
                        codon = ''.join([comp[nuc] for nuc in codon])
                        return(codon)

    def getAlFreq(self, vcf_rec):
        '''Calculate allele frequency based on allelic depth'''
        self.logger.debug('Calulating allele frequency')
        try:
            if 'DP4' in vcf_rec.INFO:
                self.logger.debug('Variant with DP4 notation')
                alt = sum(vcf_rec.INFO['DP4'][2:])
                total = sum(vcf_rec.INFO['DP4'])
                alfreq = alt/float(total)
            else:
                sample = list(vcf_rec.Samples.keys())[0]
                self.logger.debug('Variant with AD notation')
                alt = vcf_rec.Samples[sample]['AD'][1]
                total = sum(vcf_rec.Samples[sample]['AD'])
                alfreq = alt/float(total)
        # Fix added for cases where AD is missing
        # Instance found when GT == ./.
        except (KeyError, ZeroDivisionError):
            self.logger.debug('Variant laccking allele split up')
            if 'AF' in vcf_rec.INFO:
                alfreq = vcf_rec.INFO['AF'][0]
            else:
                alfreq = None
        return(alfreq)

    def createReadTable(self, readset, variant):
        readtable = list()
        for reads in readset:
            #print(reads.query_name, reads.reference_start, reads.reference_end)
            ref_start = reads.reference_start
            ref_end =  reads.reference_end
            if ref_end == None or ref_start == None:
                continue
            else:
                overlap = ref_end - ref_start + 1
                center = (ref_start+ref_end)/2
                skew = (variant.POS)/ center
                if skew >= 1 :
                    divergence = skew - 1
                else:
                    divergence = 1 - skew
                readtable.append([ref_start, ref_end, overlap, divergence, variant.POS, center])
        readtable = pd.DataFrame(readtable, columns=['Start', 'Stop', 'Overlap', 'Divergence', 'VarPos', 'Center'])
        return(readtable)

    def calculateVQ(self, bam_path, vcf_rec, cut_off=0.5):
        '''Given BAM file and vcf_rec, scan the alignments to calculate variant quality, which
        is defined as mean_centrality of var within reads, times mean_overlap between reads, times
        the depth at the loci'''

        bamreader = pysam.AlignmentFile(bam_path, 'rb')
        read_set = bamreader.fetch(region='{0}:{1}-{2}'.format(vcf_rec.CHROM, vcf_rec.POS, vcf_rec.POS+1))
        read_table = self.createReadTable(read_set, vcf_rec)
        count = 0
        min_overlap = list()
        index = [var for var in range(0, len(read_table.index))]
        read_table.sort_values(by=['Divergence','Overlap'], ascending=[True, False], inplace=True)
        read_table_sorted = read_table.reindex()
        max_index, max_overlap =  next(read_table.iterrows())
        for index, values in read_table.iterrows():
            local_start = 0
            local_end = 0
            if values.Start >= max_overlap.Start:
                local_start = values.Start
            else:
                local_start = max_overlap.Start
            if values.Stop <= max_overlap.Stop:
                local_end = values.Stop
            else:
                local_end = max_overlap.Stop

            local_overlap = local_end - local_start
            if local_overlap <= cut_off * max_overlap.Overlap:
                count += 1
            min_overlap.append(local_overlap)
        min_series = pd.Series(min_overlap)
        read_table['MinOverlap'] = min_series.values
        read_table['Centrality'] = 1 - read_table['Divergence']
        mean_centrality = np.mean(read_table.Centrality)
        mean_overlap = np.mean(read_table.MinOverlap)
        depth = len(read_table.index)
        #vardepth = vcf_rec.INFO['DP'][0]
        variant_quality = mean_centrality * mean_overlap * depth
        log_vq = np.log10(variant_quality)
        return(log_vq)

    def getAnnotation(self, bed_path, vcf_path, fasta_path, out_path, bam_path):
        #Open bed reader
        bed = Bed(bed_path)
        bed_reader = bed.getExonTable()
        #Open vcf reader
        vcf = Reader(vcf_path, fasta_path)
        vcf.readheader()
        vcf_reader = vcf.readvcf()
        #Set output files
        out_path = os.path.abspath(out_path)
        file_name = os.path.splitext(os.path.basename(vcf_path))[0]
        out_file = '{0}/{1}_annotated.vcf'.format(out_path, file_name)
        #Open vcf writer
        vcf_writer = Writer(out_file)
        try:
            vcf_rec = next(vcf_reader)
            bed_rec = next(bed_reader)
        except StopIteration:
            vcf_writer.writeHeaders(vcf.header)
            vcf_writer.closeWriter()
            return(out_file)
        #mod_temp = tempfile.NamedTemporaryFile(delete=False)
        mod_fasta =  self.getModFasta(fasta_path, vcf_path)
        mod_file = '{0}/mod_fasta.fa'.format(out_path)
        mod_out = open(mod_file, 'w')
        for sequences in mod_fasta:
            mod_out.write('>{0}\n{1}\n'.format(sequences.header,
                                                sequences.seq))
        mod_out.close()
        cds_pos = 0
        codon_range = None
        codon_change = [None, None, None]

        #Add annotation headers
        cdspos_header = self.addInfoHeader('CDSPos', type='Integer',
            description='CDS location of variant', number=1)
        vcf.header['info'].append(cdspos_header)
        vcf.info_dict['CDSPos'] = [1, 'Integer']
        ref_codon_header = self.addInfoHeader('RefCodon',
            description='Reference codon for the variant', number=1)
        vcf.header['info'].append(ref_codon_header)
        vcf.info_dict['RefCodon'] = [1, 'String']
        aapos_header = self.addInfoHeader('AAPos', type='Integer',
            description='Amino acid position altered by variant', number=1)
        vcf.header['info'].append(aapos_header)
        vcf.info_dict['AAPos'] = [1, 'Integer']
        alt_codon_header = self.addInfoHeader('AltCodon',
            description='Alternate codon for the variant', number=1)
        vcf.header['info'].append(alt_codon_header)
        vcf.info_dict['AltCodon'] = [1, 'String']
        ref_aa_header = self.addInfoHeader('RefAA',
            description='Reference amino acid', number=1)
        vcf.header['info'].append(ref_aa_header)
        vcf.info_dict['RefAA'] = [1, 'String']
        alt_aa_header = self.addInfoHeader('AltAA',
            description='Alternate amino acid', number=1)
        vcf.header['info'].append(alt_aa_header)
        vcf.info_dict['AltAA'] = [1, 'String']
        freq_header = self.addInfoHeader('Freq', type='Float',
            description='Variant allele frequency', number=1)
        vcf.header['info'].append(freq_header)
        vcf.info_dict['Freq'] = [1, 'Float']
        vqual_header = self.addInfoHeader('VFD', type='Float',
            description=('Log variant flanking depth;'
                'A measure that account for position of the base in the reads, '
                'the coverage at the location, '
                'and the average length of the overlap between the supporting reads'), number=1)
        vcf.header['info'].append(vqual_header)
        vcf.info_dict['VFD'] = [1, 'Float']
        var_header = self.addInfoHeader('Var', type='String',
            description='Variant string, with chromsome and position',
            number=1)
        vcf.header['info'].append(var_header)
        vcf.info_dict['Var'] = [1, 'String']
        loc_header = self.addInfoHeader('Exon', type='String',
            description='Exon number of the variant', number=1)
        vcf.header['info'].append(loc_header)
        vcf.info_dict['Exon'] = [1, 'String']
        gene_header = self.addInfoHeader('Gene', type='String',
            description='Gene name', number=1)
        vcf.header['info'].append(gene_header)
        vcf.info_dict['Gene'] = [1, 'String']
        #Write vcf headers
        vcf_writer.writeHeaders(vcf.header)
        #print('Annotating : {0}'.format(vcf_path))
        vcf_count = 0
        while True:
            try:
                if vcf_rec.UID in range(bed_rec.uidStart, bed_rec.uidStop):
                    #print('Annotating vcf record in bed region')
        #Strand addition
                    if bed_rec.strand == '+':
                        var_cds = cds_pos + vcf_rec.UID - bed_rec.uidStart + 1
                    else:
                        var_cds = bed_rec.uidStop - vcf_rec.UID - cds_pos 
                    triplet, aa_pos, codon_pos, exon, gene, strand =\
                                                    self.getTriplet(bed,
                                                    fasta_path, vcf_rec.UID,
                                                    var_cds)
                    alt_triplet = self.getAltTriplet(bed, mod_file, vcf_rec,
                                                    codon_pos, strand)
                    ref_aa = self.getAA(triplet)
                    if alt_triplet == None:
                        print(vcf_rec.CHROM, vcf_rec.POS, vcf_rec.REF[0], vcf_rec.ALT[0], vcf_rec.Samples)
                        print(triplet, aa_pos, codon_pos, exon, gene, strand)
                    alt_aa = self.getAA(alt_triplet)
                    allele_freq = self.getAlFreq(vcf_rec)
                    log_vq = self.calculateVQ(bam_path, vcf_rec)
                    if len(vcf_rec.REF[0]) > 1 or len(vcf_rec.ALT[0]) > 1:
                        #Add codon pos to info
                        vcf_rec.INFO['CDSPos'] = [None]
                        #Add triplet to info
                        vcf_rec.INFO['RefCodon'] = [None]
                        #Add Amino acid position
                        vcf_rec.INFO['AAPos'] = [None]
                        #Add Alt codon Amino acid position
                        vcf_rec.INFO['AltCodon'] = [None]
                        #Add Ref Amino acid
                        vcf_rec.INFO['RefAA'] = [None]
                        #Add Alt Amino acid
                        vcf_rec.INFO['AltAA'] = [None]
                        #Add Allele frequency
                        vcf_rec.INFO['Freq'] = [allele_freq]
                        #Add Log variant flanking depth
                        vcf_rec.INFO['VFD'] = [log_vq]
                        #Add Exonic location
                        vcf_rec.INFO['Exon'] = [None]
                        vcf_rec.INFO['Exon'] = [None]
                        #Add Gene
                        vcf_rec.INFO['Gene'] = [None]
                        #Add variants
                        vcf_rec.INFO['Var'] = ['{0}:{1}{2}{3}'.format(
                                            vcf_rec.CHROM, vcf_rec.REF[0],
                                            vcf_rec.POS, vcf_rec.ALT[0])]
                    else:
                        #Add codon pos to info
                        vcf_rec.INFO['CDSPos'] = [var_cds]
                        #Add triplet to info
                        vcf_rec.INFO['RefCodon'] = [triplet]
                        #Add Amino acid position
                        vcf_rec.INFO['AAPos'] = [int(aa_pos)]
                        #Add Alt codon Amino acid position
                        vcf_rec.INFO['AltCodon'] = [alt_triplet]
                        #Add Ref Amino acid
                        vcf_rec.INFO['RefAA'] = [ref_aa]
                        #Add Alt Amino acid
                        vcf_rec.INFO['AltAA'] = [alt_aa]
                        #Add Allele frequency
                        vcf_rec.INFO['Freq'] = [allele_freq]
                        #Add log variant flanking depth
                        vcf_rec.INFO['VFD'] = [log_vq]
                        #Add Exonic location
                        vcf_rec.INFO['Exon'] = [exon]
                        #Add Gene
                        vcf_rec.INFO['Gene'] = [gene]
                        #Add Variant
                        vcf_rec.INFO['Var'] = ['{0}:{1}{2}{3}'.format(
                                            gene, ref_aa, int(aa_pos), alt_aa)]

                    vcf_writer.writeRecords(vcf_rec)
                    vcf_count += 1
                    #print('Vcf record {0} written'.format(vcf_count))
                    #print(vcf_rec.CHROM, vcf_rec.POS, vcf_rec.REF, vcf_rec.ALT)

                    vcf_rec = next(vcf_reader)
                    #print('Next records are :')
                    #print(' VCF : {0} {1} {2} {3} {4}'.format(vcf_rec.CHROM, vcf_rec.POS, vcf_rec.REF, vcf_rec.ALT, vcf_rec.UID))
                    #print(' BED : {0} {1} {2} {3}'.format(bed_rec.chrom, bed_rec.gene, bed_rec.uidStart, bed_rec.uidStop))

                elif vcf_rec.UID >= bed_rec.uidStop:
                    #print('Skipping bed record')
                    #print(bed_rec.chrom, bed_rec.uidStart, bed_rec.uidStop)
                    #print('Vcf at position : {0}'.format(vcf_rec.UID))
                    current_bed_rec = bed_rec.uidStart
                    current_bed_chrom = bed_rec.chrom
                    cds_pos += bed_rec.length + 1
                    bed_rec = next(bed_reader)
                    #print('Next records are :')
                    #print(' VCF : {0} {1} {2} {3} {4}'.format(vcf_rec.CHROM, vcf_rec.POS, vcf_rec.REF, vcf_rec.ALT, vcf_rec.UID))
                    #print(' BED : {0} {1} {2} {3}'.format(bed_rec.chrom, bed_rec.gene, bed_rec.uidStart, bed_rec.uidStop))

                    # Fix added to account for MT chromsome and gene issue
                    # Codon pos was not being reset when MT gene changed
                    # TODO: Rethink strategy to handle such issues
                    # This needs to change when doing whole genome analysis
                    # because all whole genome bacterial instances will be
                    # similar to MT situation
                    if (current_bed_rec > bed.uids[bed_rec.chrom]
                        and bed_rec.chrom != bed_rec.gene):
                        cds_pos = 0
                    # If magnitude of bed uid changed by an order,
                    # it indicates a change of chromosome, hence reset
                    # CDS pos
                    if current_bed_rec < bed.uids[bed_rec.chrom]:
                        cds_pos = 0
                elif vcf_rec.UID < bed_rec.uidStart:
                    #print('Annotating vcf record in intronic region')
                    allele_freq = self.getAlFreq(vcf_rec)
                    log_vq = self.calculateVQ(bam_path, vcf_rec)
                    #Add codon pos to info
                    vcf_rec.INFO['CDSPos'] = [None]
                    #Add triplet to info
                    vcf_rec.INFO['RefCodon'] = [None]
                    #Add Amino acid position
                    vcf_rec.INFO['AAPos'] = [None]
                    #Add Alt codon Amino acid position
                    vcf_rec.INFO['AltCodon'] = [None]
                    #Add Ref Amino acid
                    vcf_rec.INFO['RefAA'] = [None]
                    #Add Alt Amino acid
                    vcf_rec.INFO['AltAA'] = [None]
                    #Add Allele frequency
                    vcf_rec.INFO['Freq'] = [allele_freq]
                    #Add log variant flanking depth
                    vcf_rec.INFO['VFD'] = [log_vq]
                    #Add Exonic location
                    vcf_rec.INFO['Exon'] = ['Intron']
                    #Add Gene
                    vcf_rec.INFO['Gene'] = [None]
                    #Add Variant
                    vcf_rec.INFO['Var'] = ['{0}:{1}{2}{3}'.format(
                                        vcf_rec.CHROM, vcf_rec.REF[0],
                                        int(vcf_rec.POS), vcf_rec.ALT[0])]
                    vcf_writer.writeRecords(vcf_rec)
                    vcf_count += 1
                    #print('Vcf record {0} written'.format(vcf_count))
                    #yield vcf_rec
                    vcf_rec = next(vcf_reader)
                    #print('Next records are :')
                    #print(' VCF : {0} {1} {2} {3} {4}'.format(vcf_rec.CHROM, vcf_rec.POS, vcf_rec.REF, vcf_rec.ALT, vcf_rec.UID))
                    #print(' BED : {0} {1} {2} {3}'.format(bed_rec.chrom, bed_rec.gene, bed_rec.uidStart, bed_rec.uidStop))

            except StopIteration:
                break
        return(out_file)
