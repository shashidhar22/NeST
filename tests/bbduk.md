# Trimming Fastq files using QualCheck

Listed below are the test cases for the different modules in the QualCheck class

## Creating an object of the QualCheck class:

The class is initialized with the paths to the following input parameters:
1. bbduk_path (str) - Path to the BBDuk executable
2. adp_path (str) - Path to the fasta file with adapter sequences
3. out_path (str) - Path to the output directory
4. java (str) - Java executable path

Create an object:

```{python}
>>> from pyamd.bbduk import QualCheck
>>> trimmer = QualCheck('bbduk.sh', 'tests/bbduk/ada.fa', 'tests/bbduk', 'java')
```

## Running QualCheck on Fastq files

Given a Fastq files, the bbduk module can be used to trim reads of low quality
remove adapter contamination. The module return the path to the quality trimmed
Fastq paths, and the returncode for the process execution. The input parameters
are:
1. rone_path (str) - Path to R1 fastq file
2. rtwo_path (str) - Path to R2 fastq file


```{python}
>>> orone_path, ortwo_path, returncode = trimmer.bbduk(
  'tests/bbduk/SRR6463555_1.fastq', 'tests/bbduk/SRR6463555_2.fastq')
```
