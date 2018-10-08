import os 
import logging
import subprocess
from nest.prepinputs import Prepper
from pprint import pprint


#Creating logger for nest
logger = logging.getLogger('NeST TB Study Prep')
logger.setLevel(logging.DEBUG)
# Creating a console handler to log info messages
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
# Create formatter and add it to the handlers
formatter = logging.Formatter('{asctime} - {name} - {levelname} - {message}', style="{")
ch.setFormatter(formatter)
# Add the handlers to the logger
logger.addHandler(ch)
#Create file and console handlers for MaRS
logger.info('Gathering input information from input path. And dowloading sample from Bioproject PRJNA271805')
prepper_sra = Prepper('fq/ColmanEtAl/SRR_Acc_List.txt', 'fastq-dump')
config = prepper_sra.prepInputs()
logger.info('{0} SRA files from PRJNA271805 downloading. Merging replicates'.format(len(config)))
sra_run = open('fq/ColmanEtAl/SraRunTable.txt')
sra_run_dict = dict()
for lineno, lines in enumerate(sra_run):
    if lineno == 0:
        continue
    lines = lines.strip().split('\t')
    sra_run_dict[lines[5]] = lines[7]

sra_file_dict = dict()
for samples in config:
    if '_1.fastq.gz' in config[samples].files[0]:
        r1 = config[samples].files[0]
        r2 = config[samples].files[1]
    else:
        r1 = config[samples].files[1]
        r2 = config[samples].files[0]
    name = sra_run_dict[samples]
    try:
        sra_file_dict[name][0].append(r1)
        sra_file_dict[name][1].append(r2)
    except KeyError:
        sra_file_dict[name] = [[r1], [r2]]

for samples in sra_file_dict:
    r1s = sra_file_dict[samples][0]
    r2s = sra_file_dict[samples][1]
    logger.info('Merging {0} into {1}'.format(';'.join(r1s), samples))
    out_r1 = 'fq/ColmanEtAl/{0}_1.fq.gz'.format(samples)
    out_r2 = 'fq/ColmanEtAl/{0}_2.fq.gz'.format(samples)
    z1 = ['cat'] + r1s + ['>', out_r1]
    z2 = ['cat'] + r2s + ['>', out_r2]
    z1r = subprocess.Popen(' '.join(z1), shell=True)
    z1r.wait()
    z2r = subprocess.Popen(' '.join(z2), shell=True)
    z2r.wait()
    allfastq = r1s + r2s
    for files in allfastq:
        os.remove(files)
