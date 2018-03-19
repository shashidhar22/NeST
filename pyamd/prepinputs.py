import os
import sys
import re
import glob
import logging
import subprocess
from collections import namedtuple
from pyamd.parsers.fastq import Fastq
from itertools import groupby


class Identifier:

    def __init__(self, record):
        self.rec = record


    def isIlluminaOld(self):
        #@HWUSI-EAS100R:6:73:941:1973#0/1
        header_regex = re.compile('@\w+-?\w+:\d+:\d+:\d+:\d+#\d*')
        match = re.fullmatch(header_regex, self.rec.header)
        if match != None:
            return(True)
        else:
            return(False)


    def isIlluminaNew(self):
        #@D00468:24:H8ELMADXX:1:1101:1470:2237 1:N:0:2
        header_regex = re.compile('@\w+-?\w+:\d+:\w+-?\w+:\d+:\d+:\d+:\d+\s\d:\w+:\w+:\w*')
        match = re.fullmatch(header_regex, self.rec.header)
        if match != None:
            return(True)
        else:
            return(False)


    def isSraOld(self):
        #@SRR037455.1 HWI-E4_6_30ACL:4:1:0:29 length=35
        #@SRR902931.1 HWI-ST1384:61:D1DJ4ACXX:8:1101:1240:2015 length=50
        header_regex = re.compile('@\w+\.?\w? \w+-?\w+:\d+:\d+:\d+:\d+ length=\d+')
        match = re.fullmatch(header_regex, self.rec.header)
        if match != None:
            return(True)
        else:
            return(False)

    def isSraNew(self):
        #@SRR1873770.5 DH1DQQN1:437:HACT2ADXX:1:2204:8270:58140 length=150
        header_regex = re.compile('@\w+\.?\w? \w+-?\w+:\d+:\w+:\d+:\d+:\d+:\d+ length=\d+')
        match = re.fullmatch(header_regex, self.rec.header)
        if match != None:
            return(True)
        else:
            return(False)

    def isENA(self):
        # ERR161234.14 14 length=100
        header_regex = re.compile('@[\w\.]+ \d+ length=\d+')
        match = re.fullmatch(header_regex, self.rec.header)
        if match != None:
            return(True)
        else:
            return(False)

    def isPacbio(self):
        #@m160113_152755_42135_c100906712550000001823199104291667_s1_p0/15/7044_26271
        header_regex = re.compile('@\w+/\d+/\d+_\d+')
        match = re.fullmatch(header_regex, self.rec.header)
        if match != None:
            return(True)
        else:
            return(False)

class Metrics:

    def __init__(self, fastq):
        self.fastq = fastq

    def avgReadLen(self):
        fastq_reader = Fastq(self.fastq, './', 'phred33')
        total_length = 0
        total_reads = 0
        for lines in fastq_reader.read():
            total_length += len(lines.seq)
            total_reads += 1
            if total_reads >= 100:
                break

        avg_length = total_length/float(total_reads)
        return(avg_length)

class Prepper:

    def __init__(self, input_path, sra_path):
        self.input_path = os.path.abspath(input_path)
        self.sra_path = sra_path
        self.logger = logging.getLogger('Kookaburra.prepInputs')

    def downloadSRA(self):
        out_dir = os.path.dirname(self.input_path)
        sra_list = open(self.input_path)
        for accessions in sra_list:
            self.logger.debug('Downloading : {0}'.format(accessions))
            accessions = accessions.strip()
            fqd_cmd = [self.sra_path, '--gzip', '--split-3', '-O', out_dir,
                        '-A', accessions]
            fqd_run = subprocess.Popen(fqd_cmd, shell=False,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
            fqd_run.wait()
            if fqd_run.returncode != 0:
                self.logger.error('Could not download {0}'.format(accessions))
                self.logger.error(' '.join(fqd_cmd))
            else:
                self.logger.info('Downladed complete: {0}'.format(accessions))
        return(out_dir)

    def getFastqPaths(self):
        filenames = list()
        for subdir, dirname, files in os.walk(self.input_path):
            for filename in files:
                if ('.fastq' in filename or '.fastq.gz' in filename or
                    'fq' in filename or 'fq.gz' in filename):
                    filepath = subdir + os.sep + filename
                    filenames.append(filepath)
        return(filenames)

    def getReadNumbers(self, file_name):
        reader = Fastq(file_name, None, None)
        read_number = 0
        for rec in reader.read():
            read_number += 1
        return(read_number)

    def prepInputs(self):
        if os.path.isfile(self.input_path):
            self.logger.info('Found SRA accession list,'
                            'Will download files from SRA')
            self.input_path = self.downloadSRA()
        files = self.getFastqPaths()
        experiment = dict()
        for fastq in files:
            reader = Fastq(fastq, './', 'phred33')
            Sample = namedtuple('Sample', ['sample', 'libname', 'library', 'files', 'prep', 'paired'])
            rec = next(reader.read())
            identifier = Identifier(rec)
            metric = Metrics(fastq)
            isIllOld =  identifier.isIlluminaOld()
            isIllNew =  identifier.isIlluminaNew()
            isSraOld = identifier.isSraOld()
            isSraNew = identifier.isSraNew()
            isPac = identifier.isPacbio()
            isENA = identifier.isENA()
            seqType = ''
            libType = ''
            sample_regex = re.compile('_r1|_r2|_?l001|_?l002|_?l003|_?l004|_R1|_R2|_L001|_?L002|_L003|_L004|_1|_2') #|L001|L002|L003|L004')
            sample = sample_regex.split(os.path.basename(fastq))[0]
            if isIllOld:
                paired_regex = re.compile('@\w+-?\w+:\d+:\d+:\d+:\d+#\d')
                lib = re.findall(paired_regex, rec.header)[0]
                paired = False
                seqType = 'Illumina'
                if metric.avgReadLen():
                    libType = 'Short'
            elif isIllNew:
                paired_regex = re.compile('@\w+-?\w+:\d+:\w+-?\w+:\d+:\d+:\d+:\d+\s')
                lib = re.findall(paired_regex, rec.header)[0]
                paired = False
                seqType = 'Illumina'
                if metric.avgReadLen():
                    libType = 'Short'
            elif isSraOld:
                paired_regex = re.compile('@\w+\.?\w? \w+-?\w+:\d+:\w+:\d+:\d+:\d+:\d+ length=\d+')
                lib = re.findall(paired_regex, rec.header)[0]
                paired = False
                seqType = 'Illumina'
                if metric.avgReadLen():
                    libType = 'Short'
            elif isSraNew:
                paired_regex = re.compile('@\w+\.?\w? \w+-?\w+:\d+:\w+:\d+:\d+:\d+:\d+ length=\d+')
                lib = re.findall(paired_regex, rec.header)[0]
                paired = False
                seqType = 'Illumina'
                if metric.avgReadLen():
                    libType = 'Short'

            elif isENA:
                paired_regex = re.compile('@[\w\.]+ \d+ length=\d+')
                lib = re.findall(paired_regex, rec.header)[0]
                paired = False
                seqType = 'Illumina'
                if metric.avgReadLen():
                    libType = 'Short'

            elif isPac:
                lib_regex = re.compile('@\w+_\d+_\d+_\w+')
                lib = re.findall(lib_regex, rec.header)[0]
                paired = False
                seqType = 'Pacbio'
                if metric.avgReadLen():
                    libType = 'Long'
            else:
                #logger.warning('Read from {0} with header : {1} does not follow any defined fastq header format.Please correct it'.format(fastq, rec_header))
                a = None
            try:
                paired = True
                #numreads = self.getReadNumbers(experiment[sample].files[0])
                experiment[sample] = Sample(sample, lib, seqType, [experiment[sample].files[0],fastq], libType, paired)
            except (KeyError, AttributeError):
                #numreads = self.getReadNumbers(fastq)
                experiment[sample] = Sample(sample, lib, seqType, [fastq], libType, paired)
        #logger.info('A total of {0} libraries were identified from the given folder {1}'.format(len(experiment), self.input_path))
        #logger.debug('The following libraries were detected in the given folder : {0}'.format(self.input_path))
        #for sample, values in experiment.items():
        #    logger.debug('Sample : {0}; Library: {1} ; Sequence type: {2} ; Files: {3} ; Library type: {4} ; Paired: {5}'.format(
        #            values.sample, values.libname, values.library, ''.join(values.files), values.prep, values.paired))
        return(experiment)

if __name__ == '__main__':
    path = os.path.abspath(sys.argv[1])
    prepper = Prepper(path)
    experiment = prepper.prepInputs()
    rone = list()
    rtwo = list()
    for study in experiment:
        for files in experiment[study].files:
            if re.findall('.*_R1.*', files):
                rone.append(files)
            else:
                rtwo.append(files)
    print(rone)
    print(rtwo)
