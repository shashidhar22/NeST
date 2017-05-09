import os
import sys
import re
import glob
import logging
from collections import namedtuple
from pyamd.readers import Fastq

logger = logging.getLogger('Assembler')

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
        header_regex = re.compile('@\w+\.?\w? \w+-\w+:\d+:\d+:\d+:\d+ length=\d+')
        match = re.fullmatch(header_regex, self.rec.header)
        if match != None:
            return(True)
        else:
            return(False)

    def isSraNew(self):
        header_regex = re.compile('@\w+\.?\w? \w+-\w+:\d+:\w+:\d+:\d+:\d+:\d+ length=\d+')
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

    def __init__(self, input_path):
        self.input_path = os.path.abspath(input_path)

    def getFastqPaths(self):
        filenames = list()
        for subdir, dirname, files in os.walk(self.input_path):
            for filename in files:
                if '.fastq' in filename or '.fastq.gz' in filename:
                    filepath = subdir + os.sep + filename
                    print(filepath)     
                    filenames.append(filepath)
        return(filenames)

    def prepInputs(self):
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
            seqType = ''
            libType = ''
            sample_regex = re.compile('_r1|_r2|r1|r2|_l001|_l002|l001|l002|l003|l004|_R1|_R2|R1|R2|_L001|_L002|L001|L002|L003|L004')
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
                paired_regex = re.compile('@\w+\.?\w? \w+-\w+:\d+:\w+:\d+:\d+:\d+:\d+ length=\d+')
                lib = re.findall(paired_regex, rec.header)[0]
                paired = False
                seqType = 'Illumina'
                if metric.avgReadLen():
                    libType = 'Short'
            elif isSraNew:
                paired_regex = re.compile('@\w+\.?\w? \w+-\w+:\d+:\w+:\d+:\d+:\d+:\d+ length=\d+')
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
                logger.warning('Read from {0} with header : {1} does not follow any defined fastq header format.Please correct it'.format(fastq, rec_header))

            try:
                paired = True
                experiment[sample] = Sample(sample, lib, seqType, [experiment[sample].files[0], fastq], libType, paired)
            except KeyError:
                experiment[sample] = Sample(sample, lib, seqType, [fastq], libType, paired)

        logger.info("The following libraries were detected in the given folder : {0}".format(self.input_path))
        for files in experiment:
            logger.info("Library: {0} ; Sequence type: {1} ; Files: {2} ; Library type: {3} ; Paired: {4}")

        return(experiment)

if __name__ == '__main__':
    path = os.path.abspath(sys.argv[1])
    prepper = Fastq(path)
    experiment = prepper.prepInputs()
