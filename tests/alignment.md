# Aligning Fastq files using BWA Mem

Listed below are the test cases for the different modules in the BWA class

## Creating an object of the BWA class:

The class is initialized with the paths to the following input parameters:
1. bwa_path (str) - Path to the BWA executable
2. out_path (str) - Path to the output path
3. ref_path (str) - Path to reference fasta file

Create an object:

```{python}
>>> from pyamd.alignment import Bwa
>>> alinger = Bwa('bwa', 'tests/alignment', 'tests/alignment/ref.fa')
```

## Running BWA mem on Fastq files

Given a set of Fastq files, the bwamem module can be used to align the reads
against the reference, using BWA-Mem. The module return the path to SAM file and
the returncode for the execution of BWA Mem

```{python}
>>> sam_path, returncode = alinger.bwamem('tests/alignment/SRR6463555_1.fastq', 'tests/alignment/SRR6463555_2.fastq')
```

# Aligning Fastq files using Bowtie2

Listed below are the test cases for the different modules in the Bowtie2 class

## Creating an object of the Bowtie class:

The class is initialized with the paths to the following input parameters:
1. bowtie_path (str) - Path to the Bowtie2 executable
2. out_path (str) - Path to the output path
3. ref_path (str) - Path to reference fasta file

Create an object:

```{python}
>>> from pyamd.alignment import Bowtie2
>>> alinger = Bowtie2('bowtie2', 'tests/alignment', 'tests/alignment/ref.fa')
```

## Running Bowtie2 on Fastq files

Given a set of Fastq files, the bowtie module can be used to align the reads
against the reference, using Bowtie2. The module return the path to SAM file and
the returncode for the execution of Bowtie2

```{python}
>>> sam_path, returncode = alinger.bowtie('tests/alignment/SRR6463555_1.fastq', 'tests/alignment/SRR6463555_2.fastq')
```

# Aligning Fastq files using BBMap

Listed below are the test cases for the different modules in the BBMap class

## Creating an object of the BBMap class:

The class is initialized with the paths to the following input parameters:
1. bbmap_path (str) - Path to the BBMap executable
2. out_path (str) - Path to the output path
3. ref_path (str) - Path to reference fasta file

Create an object:

```{python}
>>> from pyamd.alignment import BBMap
>>> alinger = BBMap('bbmap.sh', 'tests/alignment', 'tests/alignment/ref.fa')
```

## Running BBMap on Fastq files

Given a set of Fastq files, the bbmap module can be used to align the reads
against the reference, using BBMap. The module return the path to SAM file and
the returncode for the execution of BBMap

```{python}
>>> sam_path, returncode = alinger.bbmap('tests/alignment/SRR6463555_1.fastq', 'tests/alignment/SRR6463555_2.fastq')
```


# Aligning Fastq files using SNAP

Listed below are the test cases for the different modules in the SNAP class

## Creating an object of the BBMap class:

The class is initialized with the paths to the following input parameters:
1. snap_path (str) - Path to the SNAP executable
2. out_path (str) - Path to the output path
3. ref_path (str) - Path to reference fasta file

Create an object:

```{python}
>>> from pyamd.alignment import Snap
>>> alinger = Snap('snap-aligner', 'tests/alignment', 'tests/alignment/ref.fa')
```

## Running Snap on Fastq files

Given a set of Fastq files, the snap module can be used to align the reads
against the reference, using Snap. The module return the path to SAM file and
the returncode for the execution of Snap

```{python}
>>> sam_path, returncode = alinger.snap('tests/alignment/SRR6463555_1.fastq', 'tests/alignment/SRR6463555_2.fastq')
```
