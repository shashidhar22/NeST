# Summarize test cases

Listed below are the test cases for the different modules in the Prepinputs class

## Creating an object of the prepinputs class:

The class is initialized with the paths to the following input files:
1. input_path (str) - Path to the directory containing fastq files or path
to the SRA accession list.
2. sra_path (str) - Path to the faqstq-dump executable.

Create an object:

```{python}
>>> from pyamd.prepinputs import Prepper
# To create a Prepper object from SRA accession list
>>> prepper_sra = Prepper('tests/prepinputs/SRR_Acc_List.txt', 'fastq-dump')
# To create a Prepper object from Input path folder
>>> prepper_fastq = Prepper('tests/prepinputs/', None)
```


## Downloading SRA files

The downloadSRA module can be used to download Fastq files from SRA.

Downloading SRA files:

```{python}
>>> fastq_directory = prepper_sra.downloadSRA()
```

## Get a list of fastq file path from input directory

Given the path to a directory containing Fastq files, getFastqPaths returns
a list of Fastq filenames from the folder

```{python}
>>> fastq_files = prepper_fastq.getFastqPaths()
```

## Screen files with less than a minimum number of reads

The getReadPresence module checks if a given file has a minimum number of reads.
By default, the file is checked for the presence of one read. It is worth noting
that as the threshold increase for this module will affect the speed of the
pipeline. It would be advisable to set a moderate threshold.

```{python}
>>> for files in fastq_files:
>>>    presence = prepper_fastq.getReadPresence(files)
>>>    if not presence:
>>>        print('{0} did not pass minimum threshold of reads'.format(files))
```

## Extracting study information from file names

The MaRS protocol generates Fastq files with a 20 character encoded file name
that contains the following study information:
  1. Year of sample collection : 2 characters (int)
  2. Country of sample collection : 2 characters (string)
  3. City/State/Province/Site of sample collection : 2 characters (string)
  4. Day of treatment : 2 characters (int)
  5. Type of treatment : 1 character (string [A-K])
  6. Sample ID : 4 characters (int)
  7. Genus and species : 2 characters (string)
  8. Sample type : 1 character (String (B,F,P,T,S))
  9. Markers : 3 characters (8-bit code for 8 genetic markers)
  10. Replication number : 1 character (int)

The parseMaRS module uses regular expressions to tease apart this information
from the sample file name, if present, and returns a sample record

```{python}
>>> for file in fastq_files:
>>>     print(prepper_fastq.parseMaRS())
```

## Creating a study config

Given a input directory path or sra accession list, iterate through
the list and create a sample record of reach file. Each sample record is
added to a dictionary with the sample name as the key and sample record
as the value

```{python}
>>> config = prepper_sra.prepInputs()
OR
>>> config = prepper_fastq.prepInputs()
```
