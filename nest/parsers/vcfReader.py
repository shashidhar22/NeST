import re
import os
import gzip
import math
import time
import logging
from collections import namedtuple
from collections import OrderedDict

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
    def __init__(self, vcf_path, ref_path):
        ## Edit: (10/10/19); Added fasta path as an input for initializing Reader
        ## Reason: Freebayes does not have contig field in the VCF header, Contig 
        ## field is needed to create UID for variants, which are used for merging
        self.vcf_path = os.path.abspath(vcf_path)
        self.ref_path = os.path.abspath(ref_path)
        self.reader = open(self.vcf_path)
        self.info_dict = dict()
        self.format_dict = dict()

    def parseFileFormat(self, line):
        try:
            value = line.split('<')[1].split('>')[0]
        except IndexError:
            value = line.split('=')[1]
        return value

    def parseInfo(self, line):
        # Catch Info header lines; Version, source and Description
        # are optional fields
        # ##INFO=<ID=ID,Number=number,Type=type,
        # Description="description",Source="source",
        # Version="version">
        info = line.split('<')[1].split('>')[0]
        fields = ['id', 'number', 'type', 'description', 'version',
                  'source']
        info_field = namedtuple('Info', fields)
        try:
            info_field.id = info.split('ID=', 1)[1].split(',', 1)[0]
            info_field.number = info.split('Number=', 1)[1].split(',', 1)[0]
            info_field.type = info.split('Type=', 1)[1].split(',', 1)[0]
            info_field.description = info.split('Description=')[1].split('"')[1]
        except IndexError:
            print("Info field missing mandatory information")
            return -1
        try:
            info_field.source = info.split('Source=')[1].split('"')[1]
        except IndexError:
            info_field.source = None
        try:
            info_field.version = info.split('Version=')[1].split('"')[1]
        except IndexError:
            info_field.source = None
            info_field.version = None
        self.info_dict[info_field.id] = info_field.type
        return info_field

    def parseFilter(self, line):
        # Catch Filter field headers
        # ##FILTER=<ID=ID,Description="description">
        filters = line.split('<')[1].split('>')[0]
        fields = ['id', 'description']
        filter_field = namedtuple('Filter', fields)
        try:
            filter_field.id = filters.split('ID=', 1)[1].split(',', 1)[0]
            filter_field.description = filters.split('Description=')[1].split('"')[1]
        except IndexError:
            print("Filter field missing mandatory information")
            return -1
        return filter_field

    def parseFormat(self, line):
        # Catch format field headers
        # ##FORMAT=<ID=ID,Number=number,Type=type,
        # Description="description">
        formats = line.split('<')[1].split('>')[0]
        fields = ['id', 'number', 'type', 'description']
        format_field = namedtuple('Format', fields)
        try:
            format_field.id = formats.split('ID=', 1)[1].split(',', 1)[0]
            format_field.number = formats.split('Number=', 1)[1].split(',', 1)[0]
            format_field.type = formats.split('Type=', 1)[1].split(',', 1)[0]
            format_field.description = formats.split('Description=')[1].split('"')[1]
        except IndexError:
            print("Format field missing mandatory information")
            return -1
        self.format_dict[format_field.id] = format_field.type
        return format_field

    def parseAlt(self, line):
        # Catch ALT format headers
        # ##ALT=<ID=type,Description=description>
        alt = line.split('<')[1].split('>')[0]
        fields = ['id', 'description']
        alt_field = namedtuple('Alt', fields)
        try:
            alt_field.id = alt.split('ID=', 1)[1].split(',', 1)[0]
            alt_field.description = alt.split('Description=')[1].split('"')[1]
        except IndexError:
            print("Alt field missing mandatory information")
            return -1
        return alt_field

    def parseAssembly(self, line):
        # Catch assembly information
        # ##assembly=url
        assembly = line.split('=')[1]
        return assembly

    def parseContig(self, line):
        # Store conting information; Length is mandatory, URL is
        # not required
        # ##contig=<ID=ctg1,URL=ftp://somewhere.org/assembly.fa,...>
        contig = line.split('<')[1].split('>')[0]
        fields = ['id', 'url', 'assembly', 'length']
        contig_field = namedtuple('Contig', fields)
        try:
            contig_field.id = contig.split('ID=')[1].split(',', 1)[0]
        except IndexError:
            print("Contig field missing mandatory information")
            return -1
        try: 
            contig_field.url = contig.split('URL=')[1].split(',',1)[0]
        except IndexError:
            contig_field.url = None
        try: 
            contig_field.assembly = contig.split('assembly=')[1].split(',',1)[0]
        except IndexError:
            contig_field.assembly = None
        try: 
            contig_field.length = contig.split('length=')[1].split(',',1)[0]
        except IndexError:
            contig_field.length = None
        return contig_field
   
    def parseOther(self, line):
        # Store all other header formats
        other = line.split('=', 1)
        fields = ['key', 'value']
        other_fields = namedtuple('Other', fields)
        other_fields.key = other[0]
        other_fields.value = other[1]
        return other_fields

    def parseSamples(self, line):
        # Create a list of samples and fields
        fields = line[1:].split('\t')[:9]
        samples = line[1:].split('\t')[9:]
        return (fields, samples)

    def getUid(self):
        self.uids = OrderedDict()
        try:
            for index, contigs in enumerate(self.header['contig'], 1):
                try:
                    length = int(contigs.length)
                    order = 10**(int(math.log10(length)) + 1) * index
                except TypeError:
                    order = 10**(11) * index
                self.uids[contigs.id] = order
        except KeyError:
            fasta = open(self.ref_path).read().split('>')[1:]
            for index, lines in enumerate(fasta, 1):
                ids = lines.split('\n')[0]
                contig = len(''.join(lines.split('\n')[1:]))
                order = 10**(int(math.log10(contig)) + 1) * index
                self.uids[ids] = order
        return


    # Checkpoint 6/1, 4:42 PM: createInfodict method was created to store info 
    # format validation data for record parsing. Still need to write getInfo 
    # and getSampleInfo methods to validate and parse the info in each record.
    def getInfo(self, info_data):
        info_data = info_data.split(';')
        info_dict = OrderedDict()
        for info in info_data:
            key = info.split('=')[0]
            try:
                value = info.split('=')[1]
            except IndexError:
                value = [True]
            info_type = self.info_dict[key]
            if info_type == 'Integer':
                value = list(map(int, value.split(',')))
            elif info_type == 'Float':
                value = list(map(float, value.split(',')))
            elif info_type == 'String' or info_type == 'Character':
                value = value.split(',')
            info_dict[key] = value

        return info_dict

    def getSampleInfo(self, sample_lists):
        format_fields = sample_lists[0].split(':')
        sample_dict = OrderedDict()
        for samples, data in zip(self.header['samples'], sample_lists[1:]):
            sample_fields = data.split(':')
            values = OrderedDict()
            for field, value in zip(format_fields, sample_fields):
                format_type = self.format_dict[field]
                if format_type == 'Integer':
                    value = list(map(int, value.split(',')))
                elif format_type == 'Float':
                    value = list(map(float, value.split(',')))
                elif format_type == 'String' or format_type == 'Character':
                    value  = value.split(',')
                elif format_type == 'Flag':
                    value = list(map(bool, value.split(',')))
                values[field] = value
            sample_dict[samples] = values
        return sample_dict


    def readheader(self):
        line = next(self.reader, None)
        self.header = OrderedDict()
        while line:
            line = line.strip()
            #Look for format line
            if line.split('=')[0] == '##fileformat':
                self.header['fileFormat'] = self.parseFileFormat(line)
                line = next(self.reader, None)
            # Catch Info header lines; Version, source and Description
            # are optional fields
            # ##INFO=<ID=ID,Number=number,Type=type,
            # Description="description",Source="source",
            # Version="version">
            elif line.split('=',1)[0] == '##INFO':
                try: 
                    self.header['info'].append(self.parseInfo(line))
                except KeyError:
                    self.header['info'] = [self.parseInfo(line)]
                line = next(self.reader, None)
            # Catch Filter field headers
            # ##FILTER=<ID=ID,Description="description">
            elif line.split('=',1)[0] == '##FILTER':
                try:
                    self.header['filter'].append(self.parseFilter(line))
                except KeyError:
                    self.header['filter'] = [self.parseFilter(line)]
                line = next(self.reader, None)
            # Catch format field headers
            # ##FORMAT=<ID=ID,Number=number,Type=type,
            # Description="description">
            elif line.split('=',1)[0] == '##FORMAT':
                try:
                    self.header['format'].append(self.parseFormat(line))
                except KeyError:
                    self.header['format'] = [self.parseFormat(line)]
                line = next(self.reader, None)
            # Catch ALT format headers
            # ##ALT=<ID=type,Description=description>
            elif line.split('=',1)[0] == '##ALT':
                try:
                    self.header['alt'].append(self.parseAlt(line))
                except KeyError:
                    self.header['alt'] = [self.parseAlt(line)]
                line = next(self.reader, None)
            # Catch assembly information
            # ##assembly=url
            elif line.split('=',1)[0] == '##assembly':
                try:
                    self.header['assembly'].append(self.parseAssembly(line))
                except KeyError:
                    self.header['assembly'] = [self.parseAssembly(line)]
                line = next(self.reader, None)
            # Store conting information; Length is mandatory, URL is
            # not required
            # ##contig=<ID=ctg1,URL=ftp://somewhere.org/assembly.fa,...>
            elif line.split('=',1)[0] == '##contig':
                try:
                    self.header['contig'].append(self.parseContig(line))
                except KeyError:
                    self.header['contig'] = [self.parseContig(line)]
                line = next(self.reader, None)
            # Store all other information as random headers
            elif line[:2] == '##':
                try:
                    self.header['other'].append(self.parseOther(line))
                except KeyError:
                    self.header['other'] = [self.parseOther(line)]
                line = next(self.reader, None)
            # Store all field names and sample names
            elif line[:6] == '#CHROM':
                self.header['fields'], self.header['samples'] = self.parseSamples(line)
                break
            # Now parse all lines
        

    def readvcf(self):
        line = next(self.reader, None)
        while line:
            line = line.strip()
            if line[0] == '#':
                continue
            else:
                self.getUid()
                record = namedtuple('Record', self.header['fields'][:-1] + ['Samples'])
                fields = line.split('\t')
                record.CHROM = fields[0]
                record.POS = int(fields[1])
                record.UID = self.uids[fields[0]] + int(fields[1])
                record.ID = fields[2]
                record.REF = fields[3].split(',')
                record.ALT = fields[4].split(',')
                record.QUAL = float(fields[5])
                record.FILTER = fields[6]
                info = self.getInfo(fields[7])
                record.INFO = info
                sample = self.getSampleInfo(fields[8:])
                record.Samples = sample
                yield record 
                line = next(self.reader, None)
        return
