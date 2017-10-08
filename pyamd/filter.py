import os
import re
import vcf
import vcf.utils
import csv
import sys

def has_next(iterator):
    try:
        value = next(iterator)
        return(value)
    except StopIteration:
        return(None)

def exclude(var):
    interval = {'PfK13': [[50, 2100]], 'PfCRT' : [[45, 3220]], 'PfMDR1' : [[2025, 4225],
                [62, 1904]], 'PfDHPS' : [[1, 2417]], 'PfDHFR' : [[1, 1827]],
                'MT' : [[35,5945]]}
    try:
        if var.POS in range(interval[var.CHROM][0][0], interval[var.CHROM][0][1]):
            return(False)
        elif var.CHROM == 'PfMDR1' and var.POS in range(interval[var.CHROM][1][0], interval[var.CHROM][1][1]):
            return(False)
        else:
            return(True)
    except KeyError:
        return(True)

def filterer(gatk, samtools, sam_name, out_path):
    gat_reader = vcf.Reader(filename=gatk)
    sam_reader = vcf.Reader(filename=samtools)
    gat_contig = list(gat_reader.contigs.keys())
    sam_contig = list(sam_reader.contigs.keys())
    out_vcf = '{0}/{1}_variants_merged.vcf'.format(out_path, sam_name)
    merged = vcf.Writer(open(out_vcf, 'w'), sam_reader)
    if gat_contig != sam_contig:
        print('VCF does not have the same Reference')
        return
    samt = has_next(sam_reader)
    gatk = has_next(gat_reader)
    while samt != None and gatk != None:
        if sam_contig.index(samt.CHROM) == sam_contig.index(gatk.CHROM):
            if samt.POS == gatk.POS:
                if gatk.is_indel or samt.is_indel:
                    samt = has_next(sam_reader)
                    gatk = has_next(gat_reader)
                elif (len(gatk.REF) == 2 or len(gatk.ALT) == 2) or (len(samt.REF) == 2 or len(samt.ALT) == 2):
                    samt = has_next(sam_reader)
                    gatk = has_next(gat_reader)
                elif (gatk.REF and samt.REF) and (gatk.alleles and samt.alleles):
                    samt.add_info('Found', 2)
                    if not exclude(samt):
                        merged.write_record(samt)
                    samt = has_next(sam_reader)
                    gatk = has_next(gat_reader)
            elif samt.POS < gatk.POS:
                if samt.is_indel:
                    samt = has_next(sam_reader)
                elif len(samt.REF) ==  2 or len(samt.ALT) == 2:
                    samt = has_next(sam_reader)
                else:
                    samt.add_info('Found',1)
                    if not exclude(samt):
                        merged.write_record(samt)
                    samt = has_next(sam_reader)
            elif samt.POS > gatk.POS:
                if gatk.is_indel:
                    gatk = has_next(gat_reader)
                elif len(gatk.REF) ==  2 or len(gatk.ALT) == 2:
                    gatk = has_next(gat_reader)
                else:
                    gatk.add_info('Found', 1)
                    if not exclude(gatk):
                        merged.write_record(gatk)
                    gatk = has_next(gat_reader)

            #They are on the same chromosome, merge internally
        elif sam_contig.index(samt.CHROM) < sam_contig.index(gatk.CHROM):
            if samt.is_indel:
                samt = has_next(sam_reader)
            elif len(samt.REF) ==  2 or len(samt.ALT) == 2:
                samt = has_next(sam_reader)
            else:
                samt.add_info('Found',1)
                if not exclude(samt):
                    merged.write_record(samt)
                samt = has_next(sam_reader)

        elif sam_contig.index(samt.CHROM) > sam_contig.index(gatk.CHROM):
            if gatk.is_indel:
                gatk = has_next(gat_reader)
            elif len(gatk.REF) ==  2 or len(gatk.ALT) == 2:
                gatk = has_next(gat_reader)
            else:
                gatk.add_info('Found', 1)
                if not exclude(gatk):
                    merged.write_record(gatk)
                gatk = has_next(gat_reader)

    while samt != None:
        if samt.is_indel:
            samt = has_next(sam_reader)
        elif len(samt.REF) ==  2 or len(samt.ALT) == 2:
            samt = has_next(sam_reader)
        else:
            samt.add_info('Found',1)
            if not exclude(samt):
                merged.write_record(samt)
            samt = has_next(sam_reader)

    while gatk != None:
        if gatk.is_indel:
            gatk = has_next(gat_reader)
        elif len(gatk.REF) ==  2 or len(gatk.ALT) == 2:
            gatk = has_next(gat_reader)
        else:
            gatk.add_info('Found', 1)
            if not exclude(gatk):
                merged.write_record(gatk)
            gatk = has_next(gat_reader)

    return(out_vcf)

if __name__ == '__main__':
    gatk = sys.argv[1]
    samt = sys.argv[2]

    filterer(gatk, samt, 'test', './')
