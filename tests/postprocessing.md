# Samfile post processing using Samtools

Listed below are the test casese for the different modules in the Samtools class

## Creating an object of the Samtools class:

The class is initialized with the paths to the following input parameters:
1. sam_path (str) - Path to the Samtools executable path
2. bft_path (str) - Path to the Bcftools executable path
3. out_path (str) - Path to the output directory


Create an object:

```{python}
>>> from pyamd.samtools import Samtools
>>> samtools = Samtools('samtools', 'bcftools', 'tests/postprocessing')
```

## Fix mate information in SAM file using Samtools

The samtools fixmate module can be used to fix mate information in aligned
paired reads. The fixmate module creates a wrapper around the execution of the
tool and returns the path to the fixed BAM file along with the returncode for
the process execution. The input parameters for the fixmate are:
1. sam_path (str) - Path to the sample SAM file

```{python}
>>> bam_path, returncode = samtools.fixmate('tests/postprocessing/alignments/output.sam')
```

## Sort a BAM file by coordinate using Samtools

The samtools sort module can be used to sort a BAM file by reference coordinate.
The sort module creates a wrapper around the execution of the tool and returns
the path to the sorted BAM file along with the returncode for the process
execution. The input parameters for sort are:
1. bam_path (str) - Path to the sample BAM file with mate information fixed

```{python}
>>> bam_path, returncode = samtools.sort('tests/postprocessing/alignments/output_FM.bam')
```

## Remove PCR duplicate alignments from a BAM using Samtools

The samtools rmdup module can be used to remove PCR duplicates from a BAM file.
The dedup module creates a wrapper around the execution of the tool and returns
the path to the de-duplicated BAM file along with the returncode for the process
execution. The input parameters for dedup are:
1. bam_path (str) - Path to the coordinate sorted BAM file

```{python}
>>> bam_path, returncode = samtools.dedup('tests/postprocessing/alignments/output_FM_SR.bam')
```

# Samfile post processing using Picard

Listed below are the test cases for the different modules in the Picard class

## Creating an object of the QualCheck class:

The class is initialized with the paths to the following input parameters:
1. java (str) - Java executable path
2. pic_path (str) - Path to the Picard executable
3. out_path (str) - Path to the output directory


Create an object:

```{python}
>>> from pyamd.gatk import Picard
>>> picard = Picard('java', 'picard', 'tests/postprocessing')
```

## Adding read group information to SAM files using Picard

Given a SAM file, the picard module can be used to add read group information
to the BAM file. The module return the path to the read group annotated BAM
file, and the returncode for the process execution. The input parameters for
picard are:
1. bam_path (str) - Path to the PCR de-duplicated BAM file
2. sam_name (str) - Sample name

```{python}
>>> bam_path, returncode = picard.picard('tests/postprocessing/alignments/output_FM_SR_DD_.bam', 'test')
```
