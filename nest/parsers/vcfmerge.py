import os
import copy
import logging 
import itertools
from collections import namedtuple
from collections import OrderedDict
from nest.parsers.vcfReader import Reader 
from nest.parsers.vcfwriter import Writer

class Merge:

    def __init__(self, tmp_dir, vcf_dict, ref_path):
         ## Edit: (10/10/19): Added fasta path as requirement for initialization of merge
         ## Reason: Reader was modified to require fasta path, see vcfReader for more info
         self.out_path = tmp_dir
         self.vcf_dict = vcf_dict
         self.ref_path = ref_path

    def splitter(self, vcf_list):
        vcf_len  = len(vcf_list)
        if vcf_len == 1:
            return vcf_list 
        elif vcf_len != 1:
            vcf_left = self.splitter(vcf_list[:vcf_len//2])
            vcf_right = self.splitter(vcf_list[vcf_len//2:])
            vcf_list = self.merge(vcf_left, vcf_right)
            return vcf_list
        else:
            vcf_list = self.splitter(vcf_list[:vcf_len//2]) + self.splitter(vcf_list[vcf_len//2:])
            self.splitter(vcf_list)

    def mergeheader(self, lheader, rheader):
        merge_dict = copy.deepcopy(lheader)
        for headers in rheader:
            values = rheader[headers]
            if type(values) is list:
                if headers == 'fields' or headers == 'samples':
                    continue
                elif headers == 'assembly' and headers not in merge_dict:
                    merge_dict[headers] = values 
                elif headers == 'other':
                    idlist = [fields.key for fields in merge_dict[headers]]
                    for entry in values:
                        if entry.key not in idlist:
                            merge_dict[headers].append(entry)
                else:
                    try:
                        idlist = [fields.id for fields in merge_dict[headers]]
                        for entry in values:
                            if entry.id not in idlist:
                                merge_dict[headers].append(entry)
                    except KeyError:
                        merge_dict[headers] = values
        return(merge_dict)

    def mergerecords(self, lreader, rreader, lsource, rsource):
        lline = next(lreader, None)
        rline = next(rreader, None)
        while lline and rline:
            record = namedtuple('Record', lline._fields)
            if lline.UID == rline.UID:
                if not set(lline.REF).isdisjoint(set(rline.REF)) and not set(lline.ALT).isdisjoint(set(lline.ALT)):
                    record.CHROM = lline.CHROM 
                    record.POS = lline.POS
                    record.UID = lline.UID 
                    record.ID = lline.ID 
                    record.REF = list(set(lline.REF).union(set(rline.REF)))
                    record.ALT = list(set(lline.ALT).union(set(rline.ALT)))
                    record.QUAL = lline.QUAL 
                    record.FILTER = lline.FILTER
                    info_dict = copy.deepcopy(lline.INFO)
                    for key, value in rline.INFO.items():
                        if key in info_dict:
                            continue
                        else:
                            info_dict[key] = value 
                    if 'Sources' in lline.INFO and 'Sources' in rline.INFO:
                        info_dict['Confidence'][0] += rline.INFO['Confidence'][0]
                        info_dict['Sources'][0] += rline.INFO['Sources'][0]
                    elif 'Sources' in lline.INFO and 'Sources' not in rline.INFO:
                        info_dict['Confidence'][0] += len(rsource.split(','))
                        info_dict['Sources'][0] += ',{0}'.format(rsource)
                    elif 'Sources' not in lline.INFO and 'Sources' in rline.INFO:
                        info_dict['Confidence'] = [1 + rline.INFO['Confidence'][0]]
                        info_dict['Sources'] = ['{0},{1}'.format(lsource, ','.join(rline.INFO['Sources']))]
                    else:
                        info_dict['Confidence'] = [len(lsource.split(',')) + len(rsource.split(','))]
                        info_dict['Sources'] = ['{0},{1}'.format(lsource, rsource)]
                    record.INFO = info_dict
                    sample_dict = copy.deepcopy(lline.Samples)
                    for sample in rline.Samples:
                        if sample not in sample_dict:
                            continue
                        else:
                            for values in rline.Samples[sample]:
                                if values in sample_dict[sample]:
                                    continue
                                else:
                                    sample_dict[sample][values] = rline.Samples[sample][values]
                    record.Samples = sample_dict

                    lline = next(lreader, None)
                    rline = next(rreader, None)
                    yield record
                else:
                    record.CHROM = lline.CHROM 
                    record.POS = lline.POS
                    record.UID = lline.UID 
                    record.ID = lline.ID 
                    record.REF = lline.REF if len(lline.REF[0]) < len(rline.REF[0]) else rline.REF #list(set(lline.REF).union(set(rline.REF)))
                    record.ALT = lline.ALT if len(lline.ALT[0]) < len(rline.ALT[0]) else rline.ALT #list(set(lline.ALT).union(set(rline.ALT)))
                    record.QUAL = lline.QUAL 
                    record.FILTER = lline.FILTER
                    info_dict = copy.deepcopy(lline.INFO)
                    for key, value in rline.INFO.items():
                        if key in info_dict:
                            continue
                        else:
                            info_dict[key] = value 
                    if 'Sources' in lline.INFO and 'Sources' in rline.INFO:
                        info_dict['Confidence'][0] += rline.INFO['Confidence'][0]
                        info_dict['Sources'][0] += rline.INFO['Sources'][0]
                    elif 'Sources' in lline.INFO and 'Sources' not in rline.INFO:
                        info_dict['Confidence'][0] += len(rsource.split(','))
                        info_dict['Sources'][0] += ',{0}'.format(rsource)
                    elif 'Sources' not in lline.INFO and 'Sources' in rline.INFO:
                        info_dict['Confidence'] = [1 + rline.INFO['Confidence'][0]]
                        info_dict['Sources'] = ['{0},{1}'.format(lsource, ','.join(rline.INFO['Sources']))]
                    else:
                        info_dict['Confidence'] = [len(lsource.split(',')) + len(rsource.split(','))]
                        info_dict['Sources'] = ['{0},{1}'.format(lsource, rsource)]
                    record.INFO = info_dict
                    sample_dict = copy.deepcopy(lline.Samples)
                    for sample in rline.Samples:
                        if sample not in sample_dict:
                            continue
                        else:
                            for values in rline.Samples[sample]:
                                if values in sample_dict[sample]:
                                    continue
                                else:
                                    sample_dict[sample][values] = rline.Samples[sample][values]
                    record.Samples = sample_dict

                    lline = next(lreader, None)
                    rline = next(rreader, None)
                    yield record
            elif lline.UID < rline.UID:
                record = lline
                if 'Sources' not in lline.INFO:
                   record.INFO['Confidence'] = [1]
                   record.INFO['Sources'] = [lsource]
                else:
                    record.INFO['Confidence'] = lline.INFO['Confidence']
                    record.INFO['Sources'] = [','.join(lline.INFO['Sources'])]
                yield record
                lline = next(lreader, None)
            elif lline.UID > rline.UID:
                record = rline 
                if 'Sources' not in rline.INFO:
                    record.INFO['Confidence'] = [1]
                    record.INFO['Sources'] = [rsource]
                else:
                    record.INFO['Confidence'] = rline.INFO['Confidence']
                    record.INFO['Sources'] = [','.join(rline.INFO['Sources'])]
                yield record
                rline = next(rreader, None)
        while lline:
            record = lline
            if 'Sources' not in lline.INFO:
                record.INFO['Confidence'] = [1]
                record.INFO['Sources'] = ['{0}'.format(lsource)]
            else:
                record.INFO['Confidence'] = lline.INFO['Confidence']
                record.INFO['Sources'] = [','.join(lline.INFO['Sources'])]
            yield record
            lline = next(lreader, None)

        while rline:
            record = rline 
            if 'Sources' not in rline.INFO:
                record.INFO['Confidence'] = [1]
                record.INFO['Sources'] = [rsource]
            else:
                record.INFO['Confidence'] = rline.INFO['Confidence']
                record.INFO['Sources'] = [','.join(rline.INFO['Sources'])]
            yield record
            rline = next(rreader, None)

    def merge(self, vcf_left, vcf_right):
        lsource = self.vcf_dict[vcf_left[0]]
        rsource = self.vcf_dict[vcf_right[0]]
        rname = '_'.join(rsource.split(','))
        lname = '_'.join(lsource.split(','))
        out_file = '{0}/{1}_{2}_{3}.vcf'.format(self.out_path, 'tmp', lname, rname)
        self.vcf_dict[out_file] = '{0},{1}'.format(lsource, rsource)
        out_write = Writer(out_file)
        lvcf = Reader(vcf_left[0], self.ref_path)
        lvcf.readheader()
        lreader = lvcf.readvcf()
        rvcf = Reader(vcf_right[0], self.ref_path)
        rvcf.readheader()
        rreader = rvcf.readvcf()
        out_header = self.mergeheader(lvcf.header, rvcf.header)
        info_fields = ['id', 'number', 'type', 'description', 'version', 'source']
        conf_field = namedtuple('Info', info_fields)
        conf_field.id = 'Confidence'
        conf_field.number = '1'
        conf_field.type = 'Integer'
        conf_field.description = 'Number of variant callers that identified this variant'
        conf_field.version = None
        conf_field.source = None
        source_field = namedtuple('Info', info_fields)
        source_field.id = 'Sources'
        source_field.number = '1'
        source_field.type = 'String'
        source_field.description = 'Variant callers that called the variant'
        source_field.version = None
        source_field.source = None
        out_header['info'].append(conf_field)
        out_header['info'].append(source_field)
        out_write.writeHeaders(out_header)
        merge_records = self.mergerecords(lreader, rreader, lsource, rsource)
        for records in merge_records:
            out_write.writeRecords(records)
        out_write.closeWriter()
        return [out_file] 
            
