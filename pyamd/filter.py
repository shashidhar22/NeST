import os
import re
import vcf
import vcf.utils
import csv
import sys

def filterer(gatk, samtools, sam_name, out_path):
    gat_reader = vcf.Reader(filename=gatk)
    sam_reader = vcf.Reader(filename=samtools)
    out_vcf = '{0}/{1}_variants_merged.vcf'.format(out_path, sam_name)
    merged = vcf.Writer(open(out_vcf, 'w'), sam_reader)
    for gatk, samt in vcf.utils.walk_together(gat_reader, sam_reader):
        try:
            if gatk.is_indel or samt.is_indel:
                continue
            elif (len(gatk.REF) == 2 and len(gatk.ALT) == 2) or (len(samt.REF) == 2 and len(samt.ALT) == 2):
                continue

            if (gatk.REF and samt.REF) and (gatk.alleles and samt.alleles):
                samt.add_info('Found', 2)
                merged.write_record(samt)
        except AttributeError:
            if gatk == None:
                samt.add_info('Found', 1)
                merged.write_record(samt)
            elif samt == None:
                gatk.add_info('Found', 1)
                merged.write_record(gatk)
                continue
        except StopIteration:
            out_vcf = None
            break

    return(out_vcf)

if __name__ == '__main__':
    gatk = sys.argv[1]
    samt = sys.argv[2]
    kest = sys.argv[3]

    filter(gatk, samt, kest)
