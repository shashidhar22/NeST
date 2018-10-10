import re
import os
import sys
import csv
import gzip
import math
import time
import pysam
import logging
import tempfile
import itertools
import pathlib
import pandas as pd
import numpy as np
from pprint import pprint
from collections import namedtuple
from collections import OrderedDict
from nest.parsers.bed import Bed
from nest.parsers.fasta import Fasta

class Vcf:

    class Reader:
        ''' The Reader class, parses a VCF file and creates a record object for
        each variant record listed in the VCF file. The headers are formatted as
        per VCF 4.2 format. The headers are stored as a dictionary so as to
        allow easy addition of new headers to the reader object. Each variant
        record is stored as Record object with the following fields:
            1. CHROM : Chromosome information, stored as a STR
            2. POS : Position of the variant, stored as an INT
            3. ID : ID information, stored as a STR
            4. REF : Reference base, stored as a LIST of STR
            5. ALT : Alternate base, stored as a LIST of STR
            6. QUAL : Quality of variant call, stored as a FLOAT
            7. FILTER : Filter field informtion, stored as a STR
            8. Samples : Dictionary containing FORMAT field information for each
                     sample, the dictionary contains all the FORMAT fields which
                     can be accessed as a key value pair and an extra key
                     'sample', which contains the sample name
            9. INFO : Dictionary containing the INFO field information,
                     the fields are stored as a key value pair, and correspond
                     to the INFO fields described in the header
        '''
        def __init__(self, vcf_path):
            self.vcf_path = os.path.abspath(vcf_path)
            # Get header will return three values :
            # 1. rec_fields : List of field names of the header line of VCF
            # 2. samples : List of samples found in VCF file
            # 3. vcf_reader : Iterator paused at the first record in the VCF
            self.rec_fields, self.samples, self.vcf_reader = self.getHeaders()
            self.logger = logging.getLogger('NeST.VcfReader')
            if self.vcf_reader is None:
                self.logger.error('Emtpy VCF file; Please check file : {0}'.format(
                                                                self.vcf_path))
                #sys.exit()
            # Append annotation fields to record fields list
            self.rec_fields.append('Samples')
            self.rec_fields.append('UID')

            # Create a dictionary of info and format headers, along with number
            # and type of value they store, for quick reference.
            self.info_dict = {info.id: [info.number, info.type] \
                                    for info in self.header['info']}
            self.format_dict = {formats.id : [formats.number, formats.type] \
                                    for formats in self.header['format']}
            # Generate UID for eac contig listed in the VCF
            self.getUid()

        def getHeaders(self):
            '''Iterates through the VCF headers and groups them as VCF4.2
            standards. Creates a class accessible dictionary with each header
            group as key and list of header records as value. Each header record
            is stored as namedtuple with the key value pairing of the header
            stored as fields for the namedtuple. All tool specific headers are
            stored under the key 'random' '''
            self.header = OrderedDict()
            vcf_file = open(self.vcf_path)
            # Check to see if the file is empty. If file is emtpy return None
            try:
                rec = next(vcf_file).strip()
            except StopIteration:
                return(None, None, None)
            while True:
                try:
                    # Catch file format line form VCF 4.2 file
                    # ##fileformat=<VCF4.2>
                    if re.match('##fileformat=(?P<fileformat>.+)', rec):
                        ff = re.match('##fileformat=(?P<fileformat>.+)', rec)
                        self.header['fileFormat'] = ff.group('fileformat')
                        rec = next(vcf_file).strip()

                    # Catch Info header lines; Version, source and Description
                    # are optional fields
                    # ##INFO=<ID=ID,Number=number,Type=type,
                    # Description="description",Source="source",
                    # Version="version">
                    elif re.match(('##INFO=<ID=[\w\.]+,Number=(\d+|A|R|G|\.),'
                            'Type=\w+,Description=".+"(,Source="\w+")?(,'
                            'Version="\w+")?>'), rec):
                        info = re.match(('##INFO=<ID=(?P<ID>[\w\.]+),'
                            'Number=(?P<Number>(\d+|A|R|G|\.)),'
                            'Type=(?P<Type>\w+),'
                            'Description="(?P<Description>.+)"(,'
                            'Source="(?P<Source>\w+)")?(,'
                            'Version="(?P<Version>\w+)")?>'), rec)
                        info_field = namedtuple('Info', ['id', 'number', 'type',
                                                         'description',
                                                         'version', 'source'])
                        info_field.id = info.group('ID')
                        info_field.number = info.group('Number')
                        info_field.type = info.group('Type')
                        info_field.description = info.group('Description')
                        if info.group('Version'):
                            info_field.version = info.group('Version')
                        else:
                            info_field.version = None
                        if info.group('Source'):
                            info_field.source = info.group('Source')
                        else:
                            info_field.source = None

                        try:
                            self.header['info'].append(info_field)
                        except KeyError:
                            self.header['info'] = [info_field]
                        rec = next(vcf_file).strip()

                    # Catch Filter field headers
                    # ##FILTER=<ID=ID,Description="description">
                    elif re.match('##FILTER=<ID=\w+\.?\w+,Description=".+">',
                                    rec):

                        filters = re.match(('##FILTER=<ID=(?P<ID>\w+\.?\w+),'
                            'Description="(?P<Description>.+)">'), rec)
                        filter_field = namedtuple('Filter', ['id',
                                                             'description'])
                        filter_field.id = filters.group('ID')
                        filter_field.description = filters.group('Description')
                        try:
                            self.header['filter'].append(filter_field)
                        except KeyError:
                            self.header['filter'] = [filter_field]
                        rec = next(vcf_file).strip()

                    # Catch format field headers
                    # ##FORMAT=<ID=ID,Number=number,Type=type,
                    # Description="description">
                    elif re.match(('##FORMAT=<ID=\w+\.?\w?,'
                            'Number=(\d+|A|R|G|\.),'
                            'Type=\w+,Description=".+">'), rec):
                        formats = re.match(('##FORMAT=<ID=(?P<ID>\w+\.?\w?),'
                            'Number=(?P<Number>(\d+|A|R|G|\.)),'
                            'Type=(?P<Type>\w+),'
                            'Description="(?P<Description>.+)">'), rec)
                        format_field = namedtuple('Format', ['id', 'number',
                                                         'type',
                                                         'description'])
                        format_field.id = formats.group('ID')
                        format_field.number = formats.group('Number')
                        format_field.type = formats.group('Type')
                        format_field.description = formats.group('Description')
                        try:
                            self.header['format'].append(format_field)
                        except KeyError:
                            self.header['format'] = [format_field]
                        rec = next(vcf_file).strip()

                    # Catch ALT format headers
                    # ##ALT=<ID=type,Description=description>
                    elif re.match(('##ALT=<ID=(DEL|INS|DUP|INV|CNV|DUP:TANDEM|,'
                            'DUP:ME|INS:ME|\*),Description=".+">'), rec):
                        filters = re.match(('##ALT=<ID=(?P<ID>(\w+\[:\.]?\w?|\*)),'
                            'Description="(?P<Description>.+)">'), rec)
                        filter_field = namedtuple('Alt', ['id', 'description'])
                        filter_field.id = filters.group('ID')
                        filter_field.description = filters.group('Description')
                        try:
                            self.header['alt'].append(filter_field)
                        except KeyError:
                            self.header['alt'] = [filter_field]
                        rec = next(vcf_file).strip()

                    # Catch assembly information
                    # ##assembly=url
                    elif re.match('##assembly=.*', rec):
                        assembly = re.match('##assembly=(?P<assembly>.*)')
                        self.header['assembly'] = assembly.group('assembly')
                        rec = next(vcf_file).strip()

                    # Store conting information; Length is mandatory, URL is
                    # not required
                    # ##contig=<ID=ctg1,URL=ftp://somewhere.org/assembly.fa,...>
                    elif re.match('##contig=<ID=.+,length=\d+(,URL=.*)?>',
                                rec):
                        contig = re.match(('##contig=<ID=(?P<ID>.+),'
                            'length=(?P<length>\d+)(,URL=(?P<url>.*))?>'), rec)
                        contig_field = namedtuple('Contig', ['ID', 'length',
                                                             'url'])
                        contig_field.id = contig.group('ID')
                        contig_field.length = contig.group('length')
                        if contig.group('url'):
                            contig_field.url = contig.group('url')
                        else:
                            contig_field.url = None
                        try:
                            self.header['contig'].append(contig_field)
                        except KeyError:
                            self.header['contig'] = [contig_field]
                        rec = next(vcf_file).strip()

                    # Store all other information as random headers
                    elif re.match('##.+=.+', rec):
                        random = re.match('##(?P<field>[^=]+)=(?P<value>.+)', rec)
                        random_field = namedtuple('Random', ['field', 'value'])
                        random_field.field = random.group('field')
                        random_field.value = random.group('value')
                        try:
                            self.header['random'].append(random_field)
                        except KeyError:
                            self.header['random'] = [random_field]
                        rec = next(vcf_file).strip()

                    # Store all field names and sample names
                    elif re.match('#CHROM', rec):
                        fields = rec[1:].strip().split('\t')
                        rec_field = fields[:8]
                        samples = fields[9:]
                        #rec = next(vcf_file).strip()
                        return(rec_field, samples, vcf_file)
                except StopIteration:
                    break

            return

        def validateInfo(self, field, value, ref, alt):
            '''Validate Info header field and set values data structure as per
            VCF 4.2 guidelines. Check how many genotypes are available for a
            genotype call:
                If info number is 'A', the no. of values must be equal to no. of
                alternate alleles.
                If info number is 'R', the no.  of values must be equal to total
                no. of alleles.
                If info number is 'G', the no. of values must be equal to the
                no. of genotypes. Set to 3 for diploid analysis.
                If info number is '.', it is unbounded.
                If info number is '0', it is a FLAG field.
                If info is any other number it indicates the no. of values.'''
            if len(ref) + len(alt) == 2:
                gt = 3 #Temp fix, TODO: permanent way to determine ploidy
            elif len(ref) + len(alt) == 3:
                gt = 6
            elif len(ref) + len(alt) == 4:
                gt = 8
            info_number = self.info_dict[field][0]
            info_type = self.info_dict[field][1]
            # Split by ',' to tokenize by info field separator
            values = value.split(',')
            if info_number == 'A' and len(values) != len(alt):
                self.logger.error('Illegal INFO format for {0}'.format(field) +
                    '; Contains {0} entries'.format(len(values)) +
                    '; Should contian {0} entries'.format(len(alt)))
            elif info_number == 'R' and len(values) != len(ref) + len(alt):
                self.logger.error('Illegal INFO format for {0}'.format(field) +
                    '; Contains {0} entires'.format(len(values)) +
                    '; Should contain {0} entries'.format(len(ref)))
            elif info_number == 'G' and len(values) != gt:
                self.logger.error('Illegal INFO format for {0}'.format(field) +
                    '; Contains {0} entries'.format(len(values)) +
                    '; Should contain {0} entries'.format(gt))
            elif info_number == '.':
                self.logger.warning('Info number unbounded')
            elif info_number == '0' and (len(values) != 0 or
                                        type(values[0]) != bool):
                self.logger.warning('Illegal INFO format for {0}'.format(field) +
                    '; Contains {0} entries'.format(len(values)) +
                    '; Should contain no entries')
            elif info_number != 'R' and\
                info_number != 'A' and\
                info_number != 'G' and\
                info_number != '.' and int(info_number) != len(values):
                self.logger.info('Illegal INFO format for {0}'.format(field) +
                    '; Contains {0} entries'.format(len(values)) +
                    '; Should contain {0} entries'.format(info_number))
            else:
                self.logger.debug('Info number validated')
            # Type cast values based on header information
            try:
                # If info_number is in R,A,G,. groups, the values are stored as
                # a list of lists to maintain value consistency
                if info_number in ['R','A','G', '.']:
                    if info_type == 'Integer':
                        value = [[
                        int(val) if val != 'NaN' else None for val in\
                         value.split(',')]]
                    elif info_type == 'Float':
                        value = [[
                        float(val) if val != 'NaN' else None for val in\
                         value.split(',')]]
                    elif info_type == 'String':
                        value = [[
                        val if val != 'NaN' else None for val in\
                         value.split(',')]]
                elif int(info_number) > 1:
                    # If info_number is an integer >1, the values are stored as a
                    # list of lists to maintain value consistency
                    if info_type == 'Integer':
                        value = [[
                        int(val) if val != 'NaN' else None for val in\
                         value.split(',')]]
                    elif info_type == 'Float':
                        value = [[
                        float(val) if val != 'NaN' else None for val in\
                         value.split(',')]]
                    elif info_type == 'String':
                        value = [[
                        val if val != 'NaN' else None for val in\
                         value.split(',')]]
                # If info_number is 1, type cast and store value in a list
                elif int(info_number) == 1:
                    if info_type == 'Integer':
                        if value == 'NaN':
                            value = [None]
                        else:
                            value = [int(value)]
                    elif info_type == 'Float':
                        if value == 'NaN':
                            value = [None]
                        else:
                            value = [float(value)]
                    elif info_type == 'String':
                        if value == 'NaN':
                            value = [None]
                        else:
                            value = [value]
                # If info_number is 0, Store False
                elif int(info_number) == '0':
                    value = [True]
            except (ValueError, TypeError):
                self.logger.error('Value specified does not conform to any Info'
                            'classification')
                self.logger.error('Field : {0} ; Value : {1}'.format(field, value))
            return(value)

        def validateFormat(self, field, value, ref, alt):
            '''Validate Format header field and set values data structure as per
            VCF 4.2 guidelines.
            Check how many genotypes are available for a genotype call
            self.gt contains the number value for gt header.
            1. If info number is 'A', the no. of values must be equal to no. of
               alternate alleles.
            2. If info number is 'R', the no.  of values must be equal to total
               no. of alleles.
            3. If info number is 'G', the no. of values must be equal to the
               no. of genotypes.
            4. If info number is '.', it is unbounded
            5. If info number is '0', it is a FLAG field
            6. If info is any other number it indicates the no. of values'''

            if len(ref) + len(alt) == 2:
                gt = 3 #Temp fix, TODO: permanent way to determine ploidy
            elif len(ref) + len(alt) == 3:
                gt = 6
            elif len(ref) + len(alt) == 4:
                gt = 8
            #Identify number and type of value for format field
            format_number = self.format_dict[field][0]
            format_type = self.format_dict[field][1]
            # Split by ',' to tokenize by info field separator
            values = value.split(',')
            if format_number == 'A' and len(values) != len(alt):
                self.logger.error('Illegal FORMAT format for {0}'.format(field) +
                    '; Contains {0} entries'.format(len(values)) +
                    '; Should contian {0} entries'.format(len(alt)))
            elif format_number == 'R' and len(values) != (len(ref) +len(alt)):
                self.logger.error('Illegal FORMAT format for {0}'.format(field) +
                    '; Contains {0} entires'.format(len(values)) +
                    '; Should contain {0} entries'.format(len(ref)+len(alt)))
            elif format_number == 'G' and len(values) != gt:
                self.logger.error('Illegal FORMAT format for {0}'.format(field) +
                    '; Contains {0} entries'.format(len(values)) +
                    '; Should contain {0} entries'.format(gt))
            elif format_number == '.':
                self.logger.warning('Format number unbounded')
            elif format_number == '0' and len(values) != 0:
                self.logger.error('Illegal INFO format for {0}'.format(field) +
                    '; Contains {0} entries'.format(len(values)) +
                    '; Should contain no entries')
            elif format_number != 'R' and\
                format_number != 'A' and\
                format_number != 'G' and\
                format_number != '.' and int(format_number) != len(values):
                self.logger.error('Illegal INFO format for {0}'.format(field) +
                    '; Contains {0} entries'.format(len(values)) +
                    '; Should contain {0} entries'.format(format_number))
            else:
                self.logger.debug('Format number validated')
            try:
                # If format_number is in R,A,G,. groups, the values are stored as
                # a list of lists to maintain value consistency
                if format_number in ['R','A','G', '.']:
                    if format_type == 'Integer':
                        value = [[
                        int(val) if val != 'NaN' else None for val in\
                         value.split(',')]]
                    elif format_type == 'Float':
                        value = [[
                        float(val) if val != 'NaN' else None for val in\
                         value.split(',')]]
                    elif format_type == 'String':
                        value = [[
                        val if val != 'NaN' else None for val in\
                         value.split(',')]]
                # If format_number is an integer >1, the values are stored as a
                # list of lists to maintain value consistency
                elif int(format_number) > 1:
                    if format_type == 'Integer':
                        value = [[
                        int(val) if val != 'NaN' else None for val in\
                         value.split(',')]]
                    elif format_type == 'Float':
                        value = [[
                        float(val) if val != 'NaN' else None for val in\
                         value.split(',')]]
                    elif format_type == 'String':
                        value = [[
                        val if val != 'NaN' else None for val in\
                         value.split(',')]]
                # If format_number is 1, type cast and store value in a list
                elif int(format_number) == 1:
                    if format_type == 'Integer':
                        if value == 'NaN':
                            value = [None]
                        else:
                            value = [int(value)]
                    elif format_type == 'Float':
                        if value == 'NaN':
                            value = [None]
                        else:
                            value = [float(value)]
                    elif format_type == 'String':
                        if value == 'NaN':
                            value = [None]
                        else:
                            value = [value]
                # If format_number is 0, Store False
                elif int(format_number) == 0:
                    value = [True]
            except ValueError:
                self.logger.error('Value specified does not conform to any Format'
                            'classification')
                self.logger.error('Field : {0} ; Value : {1}'.format(field, value))
            return(value)

        def getUid(self):
            '''Create unique ID for each contig in VCF file'''
            self.uids = OrderedDict()
            for index, contigs in enumerate(self.header['contig'],1):
                length = int(contigs.length)
                order = 10**(int(math.log10(length)) + 1) * index
                self.uids[contigs.id] = order
            return

        def parseRecords(self, vcf_rec):
            '''Parse a VCF record and return a Record object, for each record in
            the VCF file'''
            # Store record values in namedtuple
            # PfCRT\t165\t.\tG\tT\t161.80\t.\t\
            # AC=2;AF=1.00;AN=2;DP=7;ExcessHet=3.0103;FS=0.000;MLEAC=2;
            # MLEAF=1.00;MQ=60.00;QD=23.11;SOR=0.941\tGT:AD:DP:GQ:PL\t\
            # 1/1:0,7:7:21:190,21,0
            record = namedtuple('Record', self.rec_fields)
            var_rec = vcf_rec
            record.CHROM = var_rec[0]
            record.POS = int(var_rec[1])
            record.UID = self.uids[var_rec[0]] + int(var_rec[1])
            record.ID = var_rec[2]
            record.REF = var_rec[3].split(',')
            record.ALT = var_rec[4].split(',')
            record.QUAL = float(var_rec[5])
            record.FILTER = var_rec[6]
            sample_info = list()
            # Get validated format values
            for samples, values in zip(self.samples, var_rec[9:]):
                fields = ['sample'] + var_rec[8].split(':')
                sample_information = OrderedDict()
                sample_information['sample'] = samples
                for field, value in zip(fields[1:], values.split(':')):
                    sample_information[field] = self.validateFormat(field,
                                                                    value,
                                                                    record.REF,
                                                                    record.ALT)
                self.logger.debug('Sample information')
                self.logger.debug(sample_information)
                sample_info.append(sample_information)
            record.Samples = sample_info
            # Get validated info values
            information = OrderedDict()
            for info in var_rec[7].split(';'):
                info = info.split('=')
                information[info[0]] = self.validateInfo(info[0], info[1],
                                                         record.REF, record.ALT)

            # Get missing inforation fields
            missing = list()
            for keys in self.header['info']:
                if keys.id in information.keys():
                    continue
                else:
                    missing.append(keys.id)
            # If there are any missing info fields, populate with null
            for fields in missing:
                value_range = self.info_dict[fields][0]
                if value_range == 'R':
                    information[fields] = [[None] * (len(record.REF) +\
                                                    len(record.ALT))]
                elif value_range == 'A':
                    information[fields] = [[None] * len(record.ALT)]
                elif value_range == 'G':
                    information[fields] = [[None] * 3]
                elif value_range == '.' or value_range == '1':
                    information[fields] = [None]
                elif int(value_range) > 1:
                    information[fields] = [[None] * int(value_range)]
                elif value_range == '0':
                    information[fields] = [False]
            self.logger.debug('Info field values')
            self.logger.debug(information)
            record.INFO = information

            return(record)

        def validateInfoFilter(self, filters):
            '''Takes a String of filters, and checks which of the filters can
            be applied as per info headers'''
            # Validate if info field filters are valid
            try:
                filters = filters.split(';')
                filter_dict = OrderedDict()
                for filter_value in filters:
                    key, value = filter_value.split('=')
                    if key in self.info_dict.keys():
                        info_type = self.info_dict[key][1]
                        if info_type == 'Integer':
                            value = int(value)
                            filter_dict[key] = value
                        elif info_type == 'Float':
                            value = float(value)
                            filter_dict[key] = value
                        elif info_type == 'Flag':
                            value == bool(value)
                            filter_dict[key] = value
                        else:
                            self.logger.error(
                                     'String field cannot be used as a filter')
                            a = None
                return(filter_dict)
            except AttributeError:
                return(None)

        def validateFormatFilter(self, filters):
            '''Takes a String of filters, and checks which of the filters can
            be applied as per info headers'''
            # Validate if format field filters are valid
            try:
                filters = filters.split(';')
                filter_dict = OrderedDict()
                for filter_value in filters:
                    key, value = filter_value.split('=')
                    if key in self.format_dict.keys():
                        format_type = self.format_dict[key][1]
                        if format_type == 'Integer':
                            value = int(value)
                            filter_dict[key] = value
                        elif format_type == 'Float':
                            value = float(value)
                            filter_dict[key] = value
                        elif format_type == 'Flag':
                            value == bool(value)
                            filter_dict[key] = value
                        else:
                            filter_dict[key] = value
                return(filter_dict)

            except AttributeError:
                return(None)

        def read(self, qual_filter=None, info_filter=None, sample_filter=None):
            '''Iterate through VCF records, filter according to user specified
            filters and yield a generator of records'''
            info_filter = self.validateInfoFilter(info_filter)
            sample_filter = self.validateFormatFilter(sample_filter)
            if self.header == None or self.rec_fields == None:
                return
            while True:
                try:
                    vcf_rec = next(self.vcf_reader).strip().split('\t')
                    rec = self.parseRecords(vcf_rec)
                    filter_info = False
                    filter_sample = False
                    filter_qual = False
                    # if there are any info qual, filters, filter records
                    if qual_filter:
                        if rec.QUAL < qual_filter:
                            filter_qual = True
                    if info_filter:
                        for key, value in info_filter.items():
                            if type(value) is int and rec.INFO[key][0] < value:
                                self.logger.info(
                                        'Skipping variant; {0}:{1}:{2}'.format(
                                        key, value, rec.INFO[key][0]))
                                filter_info = True
                            elif (type(value) is float and
                                rec.INFO[key][0] < value):
                                self.logger.info(
                                        'Skipping variant; {0}:{1}:{2}'.format(
                                        key, value, rec.INFO[key][0]))
                                filter_info = True
                            elif (type(value) is bool and\
                                rec.INFO[key] == None):
                                self.logger.info(
                                        'Skipping variant; {0}:{1}:{2}'.format(
                                        key, value, rec.INFO[key][0]))
                                filter_info = True
                            else:
                                filter_info = False

                    # if there are any sample filters, filter records
                    if sample_filter:
                        for samples in rec.Samples:
                            if samples['GT'][0] != sample_filter['GT']:
                                pprint('Skipping variant; {0}:{1}:{2}'.format(
                                    'GT', sample_filter['GT'],samples['GT'][0]))
                                filter_sample = True
                            else:
                                filter_sample = False
                    if (not filter_info and not filter_sample and
                        not filter_qual):
                        yield(rec)
                except StopIteration:
                    break

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
            vcf = Vcf.Reader(vcf_path)
            vcf_reader = vcf.read()
            fasta = Fasta(fasta_path)
            fasta_reader = fasta.read()
            for fasta_rec in fasta_reader:
                vcf = Vcf.Reader(vcf_path)
                vcf_reader = vcf.read()
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
                    alt = sum(vcf_rec.INFO['DP4'][0][2:])
                    total = sum(vcf_rec.INFO['DP4'][0])
                    alfreq = alt/float(total)
                else:
                    self.logger.debug('Variant with AD notation')
                    alt = vcf_rec.Samples[0]['AD'][0][1]
                    total = sum(vcf_rec.Samples[0]['AD'][0])
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
            read_table_sorted = read_table.reindex(labels=index)
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
            vardepth = vcf_rec.INFO['DP'][0]
            variant_quality = mean_centrality * mean_overlap * depth
            log_vq = np.log10(variant_quality)
            return(log_vq)

        def getAnnotation(self, bed_path, vcf_path, fasta_path, out_path, bam_path):
            #Open bed reader
            bed = Bed(bed_path)
            bed_reader = bed.getExonTable()
            #Open vcf reader
            vcf = Vcf.Reader(vcf_path)
            vcf_reader = vcf.read()
            #Set output files
            out_path = os.path.abspath(out_path)
            file_name = os.path.splitext(pathlib.Path(vcf_path).name)[0]
            out_file = '{0}/{1}_annotated.vcf'.format(out_path, file_name)
            #Open vcf writer
            vcf_writer = Vcf.Writer(out_file)
            try:
                vcf_rec = next(vcf_reader)
                bed_rec = next(bed_reader)
            except StopIteration:
                vcf_writer.writeHeaders(vcf.header, vcf.rec_fields, vcf.samples)
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
            vcf.info_dict['CDSPos'] = [1, 'Interger']
            ref_codon_header = self.addInfoHeader('RefCodon',
                description='Reference codon for the variant', number=1)
            vcf.header['info'].append(ref_codon_header)
            vcf.info_dict['RefCodon'] = [1, 'String']
            aapos_header = self.addInfoHeader('AAPos', type='Integer',
                description='Amino acid position altered by variant', number=1)
            vcf.header['info'].append(aapos_header)
            vcf.info_dict['AAPos'] = [1, 'Interger']
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
            vcf_writer.writeHeaders(vcf.header, vcf.rec_fields, vcf.samples)
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
                            vcf_rec.INFO['AAPos'] = [aa_pos]
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
                                                gene, ref_aa, aa_pos, alt_aa)]

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
                                            vcf_rec.POS, vcf_rec.ALT[0])]
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

    class Merge:

        def __init__(self, fone, ftwo, out_path ):
            self.reading_list = [fone, ftwo]
            self.info_fields = list()
            self.sample_fields = list()
            self.out_path = out_path
            self.info_hone = OrderedDict()
            for header in fone.header['info']:
                if header.id not in self.info_fields:
                    self.info_fields.append(header.id)
                self.info_hone[header.id] = header.number
            self.info_htwo = OrderedDict()
            for header in ftwo.header['info']:
                if header.id not in self.info_fields:
                    self.info_fields.append(header.id)
                self.info_htwo[header.id] = header.number
            self.sample_hone = OrderedDict()
            for header in fone.header['format']:
                if header.id not in self.sample_fields:
                    self.sample_fields.append(header.id)
                self.sample_hone[header.id] = header.number
            self.sample_htwo = OrderedDict()
            for header in ftwo.header['format']:
                if header.id not in self.sample_fields:
                    self.sample_fields.append(header.id)
                self.sample_htwo[header.id] = header.number

        def has_next(self, iterator):
            try:
                value = next(iterator)
                return(value)
            except StopIteration:
                return(None)

        def merge_by_id(self, field, header_one, header_two):
            rone_info = list()
            try:
                for header in header_one[field]:
                    if header.id in rone_info:
                        continue
                    else:
                        rone_info.append(header.id)
            except KeyError:
                rone_info = list()
#            try:
#                rone_info = set([header.id for header in header_one[field]])
#            except KeyError:
#                rone_info = set()

            rtwo_info = list()
            try:
                for header in header_two[field]:
                    if header.id in rtwo_info:
                        continue
                    else:
                        rtwo_info.append(header.id)
            except KeyError:
                rtwo_info = list()
#            try:
#                rtwo_info = set([header.id for header in header_two[field]])
#            except KeyError:
#                rtwo_info = set()

            common = list()
            rones = list()
            rtwos = list()
            for values in rone_info:
                if values in rone_info and values in rtwo_info and values not in common:
                    common.append(values)
                elif values in rone_info and values not in rtwo_info and values not in rones:
                    rones.append(values)

            for values in rtwo_info:
                if values in rone_info and values in rtwo_info and values not in common:
                    common.append(values)
                elif values in rtwo_info and values not in rone_info and values not in rtwos:
                    rtwos.append(values)

#            common = rone_info & rtwo_info
#            rones = rone_info - rtwo_info
#            rtwos = rtwo_info - rone_info

            try:
                rone_dict = OrderedDict()
                for header in header_one[field]:
                    rone_dict[header.id] = header

#                rone_dict = {header.id : header for header in header_one[field]}
            except KeyError:
                rone_dict = OrderedDict()

            try:
                rtwo_dict = OrderedDict()
                for header in header_two[field]:
                    rtwo_dict[header.id] = header

#                rtwo_dict = {header.id : header for header in header_two[field]}
            except KeyError:
                rtwo_dict = OrderedDict()

            for ids in common:
                header = rone_dict[ids]
                try:
                    self.merged_header[field].append(header)
                except KeyError:
                    self.merged_header[field] = [header]

            for ids in rones:
                header = rone_dict[ids]
                try:
                    self.merged_header[field].append(header)
                except KeyError:
                    self.merged_header[field] = [header]

            for ids in rtwos:
                header = rtwo_dict[ids]
                try:
                    self.merged_header[field].append(header)
                except KeyError:
                    self.merged_header[field] = [header]

        def merge_by_field(self, field, header_one, header_two):
            rone_info = list()
            try:
                for header in header_one[field]:
                    if header.field in rone_info:
                        continue
                    else:
                        rone_info.append(header.field)
            except KeyError:
                rone_info = list()
#            try:
#                rone_info = set([header.field for header in header_one[field]])
#            except KeyError:
#                rone_info = set()

            rtwo_info = list()
            try:
                for header in header_two[field]:
                    if header.field in rtwo_info:
                        continue
                    else:
                        rtwo_info.append(header.field)
            except KeyError:
                rtwo_info = list()
#            try:
#                rtwo_info = set([header.field for header in header_two[field]])
#            except KeyError:
#                rtwo_info = set()
            common = list()
            rones = list()
            rtwos = list()
            for values in rone_info:
                if values in rone_info and values in rtwo_info:
                    common.append(values)
                elif values in rone_info and values not in rtwo_info:
                    rones.append(values)

            for values in rtwo_info:
                if values in rone_info and values in rtwo_info:
                    common.append(values)
                elif values in rone_info and values not in rtwo_info:
                    rtwos.append(values)


#            common = rone_info & rtwo_info
#            rones = rone_info - rtwo_info
#            rtwos = rtwo_info - rone_info

            try:
                rone_dict = OrderedDict()
                for header in header_one[field]:
                    rone_dict[header.field] = header
                #rone_dict = {header.field : header for header in
                #            header_one[field]}
            except KeyError:
                rone_dict = OrderedDict()
            try:
                rtwo_dict = OrderedDict()
                for header in header_two[field]:
                    rtwo_dict[header.field] = header
                #rtwo_dict = {header.field : header for header in
                #            header_two[field]}
            except KeyError:
                rtwo_dict = OrderedDict()
            for ids in common:
                header = rone_dict[ids]
                try:
                    self.merged_header[field].append(header)
                except KeyError:
                    self.merged_header[field] = [header]

            for ids in rones:
                header = rone_dict[ids]
                try:
                    self.merged_header[field].append(header)
                except KeyError:
                    self.merged_header[field] = [header]

            for ids in rtwos:
                header = rtwo_dict[ids]
                try:
                    self.merged_header[field].append(header)
                except KeyError:
                    self.merged_header[field] = [header]

        def merge_headers(self, header_one, header_two):
            self.merged_header = OrderedDict()
            seen_keys = list()
            # Iterate through the first header set and merge every header gorup
            # if header group not present in second header set, write first
            # heaer set to merged dictionary

            ## File format merge
            ff = header_one['fileFormat']
            self.merged_header['fileFormat'] = ff

            ## Info merge
            ##TEMP
            self.merge_by_id('info', header_one, header_two)

            ## Filter merge
            ##TEMP
            self.merge_by_id('filter', header_one, header_two)

            ## Format merge
            ##TEMP
            self.merge_by_id('format', header_one, header_two)

            ## Alt merge
            ##TEMP
            self.merge_by_id('alt', header_one, header_two)

            ## Contig merge
            ##TEMP
            self.merge_by_id('contig', header_one, header_two)

            ## assembly merge
            if 'assembly' in header_one:
                self.merged_header['assembly'] = header_one['assembly']
            elif 'assembly' in header_two:
                self.merged_header['assembly'] = header_two['assembly']

            ## Random merge
            ##TEMP
            self.merge_by_field('random', header_one, header_two)

        def merge_record(self, record_one, record_two, header_one, header_two):
            if record_one != None and record_two != None:
                record = namedtuple('Record', record_one._fields)
                record.CHROM = record_one.CHROM
                record.POS = record_one.POS
                record.UID = record_one.UID
                record.ID = record_one.ID
                record.REF = record_one.REF
                record.ALT = record_one.ALT
                record.QUAL = record_one.QUAL
                record.FILTER = record_one.FILTER
                sample_info = list()
                for sample_one, sample_two in zip(record_one.Samples,
                                                record_two.Samples):
                    sample_information = OrderedDict()
                    for fields in self.sample_fields:
                        if fields in self.sample_hone:
                            sample_information[fields] = sample_one[fields]
                        if fields in self.sample_htwo:
                            sample_information[fields] = sample_two[fields]
                    sample_info.append(sample_information)
                record.Samples = sample_info

            elif record_one != None and record_two == None:
                record = namedtuple('Record', record_one._fields)
                record.CHROM = record_one.CHROM
                record.POS = record_one.POS
                record.UID = record_one.UID
                record.ID = record_one.ID
                record.REF = record_one.REF
                record.ALT = record_one.ALT
                record.QUAL = record_one.QUAL
                record.FILTER = record_one.FILTER
                sample_info = list()
                for sample_one in record_one.Samples:
                    sample_information = OrderedDict()
                    for fields in self.sample_fields:
                        if fields in self.sample_hone:
                            sample_information[fields] = sample_one[fields]
                        elif fields in self.sample_htwo:
                            if self.sample_htwo[fields] == 'R':
                                value = [[None] * (len(record_one.REF) +\
                                                    len(record_one.ALT))]
                            elif self.sample_htwo[fields] == 'A':
                                value = [[None] * len(record_one.ALT)]
                            elif self.sample_htwo[fields] == 'G':
                                value = [[None] * 3]
                            elif self.sample_htwo[fields] == '0':
                                value = [False]
                            elif self.sample_htwo[fields] == '.':
                                value = [None]
                            elif int(self.sample_htwo[fields]) > 1:
                                value = [[None] * int(self.sample_htwo[fields])]
                            else:
                                value = [None]
                            sample_information[fields] = value

                    sample_info.append(sample_information)
                record.Samples = sample_info

            elif record_one == None and record_two != None:
                record = namedtuple('Record', record_two._fields)
                record.CHROM = record_two.CHROM
                record.POS = record_two.POS
                record.UID = record_two.UID
                record.ID = record_two.ID
                record.REF = record_two.REF
                record.ALT = record_two.ALT
                record.QUAL = record_two.QUAL
                record.FILTER = record_two.FILTER
                sample_info = list()
                #print(record_two.CHROM, record_two.POS, record_two.REF, record_two.ALT)
                for sample_two in record_two.Samples:
                    #print(sample_two)
                    sample_information = OrderedDict()
                    for fields in self.sample_fields:
                        #print(self.sample_htwo)
                        #print(self.sample_hone)
                        if fields in self.sample_htwo:
                            #print(fields, sample_two[fields])
                            sample_information[fields] = sample_two[fields]
                        elif fields in self.sample_hone:
                            if self.sample_hone[fields] == 'R':
                                value = [[None] * (len(record_two.REF) +\
                                                    len(record_two.ALT))]
                            elif self.sample_hone[fields] == 'A':
                                value = [[None] * len(record_two.ALT)]
                            elif self.sample_hone[fields] == 'G':
                                value = [[None] * 3]
                            elif self.sample_hone[fields] == '0':
                                value = [False]
                            elif self.sample_hone[fields] == '.':
                                value = [None]
                            elif int(self.sample_hone[fields]) > 1:
                                value = [[None] * int(self.sample_hone[fields])]
                            else:
                                value = [None]
                            #print(fields, value)
                            sample_information[fields] = value
                    sample_info.append(sample_information)
                record.Samples = sample_info
            information = OrderedDict()
            for fields in self.info_fields:
                if fields in self.info_hone and record_one != None:
                    information[fields] = record_one.INFO[fields]
                elif fields in self.info_htwo and record_two != None:
                    information[fields] = record_two.INFO[fields]
                elif fields in self.info_hone and record_one == None:
                    field_number = self.info_hone[fields]
                    if field_number == 'R':
                        value = [None] * (len(record_two.REF) +\
                                            len(record_two.ALT))
                    elif field_number == 'A':
                        value = [None] * len(record_two.ALT)
                    elif field_number == 'G':
                        value = [None] * 3
                    elif field_number == '0':
                        value = [False]
                    elif field_number == '.':
                        value = [None]
                    else:
                        value = [None] * int(field_number)
                    information[fields] = value
                elif fields in self.info_htwo and record_two == None:
                    field_number = self.info_htwo[fields]
                    if field_number == 'R':
                        value = [None] * (len(record_one.REF) +\
                                            len(record_one.ALT))
                    elif field_number == 'A':
                        value = [None] * len(record_one.ALT)
                    elif field_number == 'G':
                        value = [None] * 3
                    elif field_number == '0':
                        value = [False]
                    elif field_number == '.':
                        value = [None]
                    else:
                        value = [None] * int(field_number)
                    information[fields] = value
            record.INFO = information
            return(record)

        def merge(self):
            annotate = Vcf.Annotate()
            vcf_one = self.reading_list[0]
            vcf_two = self.reading_list[1]
            self.merge_headers(vcf_one.header, vcf_two.header)
            reader_one = vcf_one.read()
            reader_two = vcf_two.read()
            record_one = next(reader_one, None)
            record_two = next(reader_two, None)
            info_headers = [header.id for header in self.merged_header['info']]
            self.merged_header['info'].append(annotate.addInfoHeader(
                                                        'Conf', type='Integer'))
            #writer objects
            sample = os.path.basename(self.out_path.strip('/'))
            out_file = '{0}/{1}_variants_merged_annotated.vcf'.format(self.out_path,
                                                                        sample)
            writer = Vcf.Writer(out_file)
            writer.writeHeaders(self.merged_header, vcf_one.rec_fields,
                                vcf_one.samples)
            while record_one != None and record_two != None:
                if  record_one.UID == record_two.UID:
                    if ('.' in record_one.REF or '.' in record_one.ALT or
                            '.' in record_two.REF or '.' in record_two.ALT):
                        record_one = next(reader_one, None)
                        record_two = next(reader_two, None)
                    elif (len(record_one.REF) > 1 or len(record_one.ALT) > 1 or
                        len(record_two.REF) > 1 or len(record_two.ALT) > 1):
                        record_one = next(reader_one, None)
                        record_two = next(reader_two, None)
                    elif (len(record_one.REF[0]) > 1 or
                        len(record_one.ALT[0]) > 1 or
                        len(record_two.REF[0]) > 1 or
                        len(record_two.ALT[0]) > 1):
                        record_one = next(reader_one, None)
                        record_two = next(reader_two, None)
                    elif ((record_one.REF and record_two.REF) and
                            (record_one.ALT and record_two.ALT)) :
                        merged_record = self.merge_record(record_one,
                            record_two, vcf_one.header, vcf_two.header)
                        merged_record.INFO['Conf'] = [2]
                        record_one = next(reader_one, None)
                        record_two = next(reader_two, None)
                        writer.writeRecords(merged_record)
                        #yield(merged_record)

                elif record_one.UID < record_two.UID:
                    if ('.' in record_one.REF or '.' in record_one.ALT):
                        record_one = next(reader_one, None)
                    elif (len(record_one.REF) > 1 or len(record_one.ALT) > 1):
                        record_one = next(reader_one, None)
                    elif (len(record_one.REF[0]) > 1 or
                        len(record_one.ALT[0]) > 1):
                        record_one = next(reader_one, None)
                    else:
                        merged_record = self.merge_record(record_one,
                            None, vcf_one.header, vcf_two.header)
                        merged_record.INFO['Conf'] = [1]
                        record_one = next(reader_one, None)
                        writer.writeRecords(merged_record)
                        #yield(merged_record)

                elif record_one.UID > record_two.UID:
                    if ('.' in record_two.REF or '.' in record_two.ALT):
                        record_two = next(reader_two, None)
                    elif (len(record_two.REF) > 1 or len(record_two.ALT) > 1):
                        record_two = next(reader_two, None)
                    elif (len(record_two.REF[0]) > 1 or
                        len(record_two.ALT[0]) > 1):
                        record_two = next(reader_two, None)
                    else:
                        merged_record = self.merge_record(None,
                            record_two, vcf_one.header, vcf_two.header)
                        merged_record.INFO['Conf'] = [1]
                        #print(merged_record.Samples)
                        record_two = next(reader_two, None)
                        writer.writeRecords(merged_record)
                        #yield(merged_record)

            while record_one != None:
                if ('.' in record_one.REF or '.' in record_one.ALT):
                    record_one = next(reader_one, None)
                elif (len(record_one.REF) > 1 or len(record_one.ALT) > 1):
                    record_one = next(reader_one, None)
                elif (len(record_one.REF[0]) > 1 or
                    len(record_one.ALT[0]) > 1):
                    record_one = next(reader_one, None)
                else:
                    merged_record = self.merge_record(record_one,
                        None, vcf_one.header, vcf_two.header)
                    merged_record.INFO['Conf'] = [1]
                    record_one = next(reader_one, None)
                    writer.writeRecords(merged_record)
                    #yield(merged_record)

            while record_two != None:
                if ('.' in record_two.REF or '.' in record_two.ALT):
                    record_two = next(reader_two, None)
                elif (len(record_two.REF) > 1 or len(record_two.ALT) > 1):
                    record_two = next(reader_two, None)
                elif (len(record_two.REF[0]) > 1 or
                    len(record_two.ALT[0]) > 1):
                    record_two = next(reader_two, None)
                else:
                    merged_record = self.merge_record(None,
                        record_two, vcf_one.header, vcf_two.header)
                    merged_record.INFO['Conf'] = [1]
                    record_two = next(reader_two, None)
                    writer.writeRecords(merged_record)
                    #yield(merged_record)
            return(out_file)

    class Writer:

        def __init__(self, out_path):
            self.out_path = out_path
            self.vcf_writer = open(out_path, 'w')

        def writeHeaders(self, headers, fields, samples):
            for header in headers:
                if header == 'fileFormat':
                    self.vcf_writer.write('##fileformat={0}\n'.format(
                                                    headers[header]))
                elif header == 'info':
                    for info in headers['info']:
                        if info.source == None:
                            self.vcf_writer.write(('##INFO=<ID={0},'.format(
                                                                        info.id)
                                 + 'Number={0},'.format(info.number)
                                 + 'Type={0},'.format(info.type)
                                 + 'Description="{0}">\n'.format(
                                                            info.description)
                                 ))
                        else:
                            self.vcf_writer.write(('##INFO=<ID={0},'.format(
                                                                        info.id)
                                 + 'Number={0},'.format(info.number)
                                 + 'Type={0},'.format(info.type)
                                 + 'Description="{0}",'.format(info.description)
                                 + 'Source="{0}",'.format(info.source)
                                 + 'Version="{0}">\n'.format(info.version)))
                elif header == 'filter':
                    for filters in headers['filter']:
                        self.vcf_writer.write(('##FILTER=<ID={0},'.format(
                                                                filters.id)
                            + 'Description="{0}">\n'.format(filters.description)
                            ))
                elif header == 'format':
                     for formats in headers['format']:
                         #print(formats.id)
                         self.vcf_writer.write(('##FORMAT=<ID={0},'.format(
                                                                formats.id)
                            + 'Number={0},'.format(formats.number)
                            + 'Type={0},'.format(formats.type)
                            + 'Description="{0}">\n'.format(formats.description)
                            ))
                elif header == 'alt':
                    for alt in headers['alt']:
                        self.vcf_writer.write(('##ALT=<ID={0},'.format(alt.id)
                            + 'Description="{0}">\n'.format(alt.description)))
                elif header == 'assembly':
                    self.vcf_writer.write('##assembly={0}\n'.format(
                                        headers[header]))
                elif header == 'contig':
                    for contig in headers['contig']:
                        if contig.url == None:
                            self.vcf_writer.write(('##contig=<ID={0},'.format(
                                                                contig.id)
                                    + 'length={0}>\n'.format(contig.length)
                                    ))
                        else:
                            self.vcf_writer.write(('##contig=<ID={0},'.format(
                                                                contig.id)
                                    + 'length={0},'.format(contig.length)
                                    + 'URL={0}>\n'.format(contig.url)))
                elif header == 'random':
                    for rands in headers['random']:
                        self.vcf_writer.write(('##{0}={1}\n'.format(
                                    rands.field, rands.value)))
            ##Write header line
            header_line = fields[:-2]
            header_line.append('FORMAT')
            for sample in samples:
                header_line.append(sample)
            self.vcf_writer.write('#{0}\n'.format('\t'.join(header_line)))

        def writeRecords(self, records):
            rec = list()
            #print(records.Samples)
            rec.append('{0}'.format(records.CHROM))
            rec.append('{0}'.format(records.POS))
            rec.append('{0}'.format(records.UID))
            rec.append(','.join(records.REF))
            alt = ','.join(records.ALT)
            rec.append(alt)
            rec.append('{0}'.format(records.QUAL))
            rec.append(records.FILTER)
            info = list()
            for key, value in records.INFO.items():
                #print(key,value)
                # If values are absent, do not write value
                if value == None or value[0] == None:
                    #info.append('{0}={1}'.format(key, 'NaN'))
                    continue
                elif type(value[0]) == list:
                    info.append('{0}={1}'.format(key,
                                ','.join([str(val) for val in value[0]])))
                elif type(value[0]) == bool and value[0]:
                    info.append('{0}'.format(key))
                elif type(value[0]) == bool and not value[0]:
                    continue
                else:
                    info.append('{0}={1}'.format(key, value[0]))
            rec.append(';'.join(info))
            all_formats = list()
            #print(records.CHROM, records.POS)
            #print(records.Samples)
            for key in records.Samples[0]:
                if key == 'sample':
                    continue
                else:
                    all_formats.append(key)
            #rec.append(':'.join(formats))
            formats = list()
            for sample in records.Samples:
                #1/1:0,41:41:99:1535,123,0
                #1/1:1:0,41:0:41:4:99:9:1535,123,
                sam_string = list()
                for keys in all_formats:
                    #print(keys, sample[keys])
                    value = sample[keys]
                    if value[0] == None:
                        continue
                    elif type(value[0]) == list and None in value[0]:
                        continue
                    elif type(value[0]) == list and len(value[0]) == 1:
                        sam_string.append(str(value[0][0]))
                        if keys not in formats:
                            formats.append(keys)
                    elif type(value[0]) == list and len(value[0]) > 1 :
                        value = ','.join([
                        str(val) if val != None else 'NaN'for val in value[0]])
                        sam_string.append(value)
                        if keys not in formats:
                            formats.append(keys)
                    else:
                        if value[0] == None:
                            continue
                        else:
                            sam_string.append(str(value[0]))
                            if keys not in formats:
                                formats.append(keys)
                #print(':'.join(formats))
                #print(':'.join(sam_string))
                rec.append(':'.join(formats))
                rec.append(':'.join(sam_string))
            self.vcf_writer.write('{0}\n'.format('\t'.join(rec)))
            return

        def closeWriter(self):
            self.vcf_writer.close()
            return
