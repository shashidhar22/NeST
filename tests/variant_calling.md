# Variant calling  using GATK

Listed below are the test cases for the different modules in the GenAnTK class

## Creating an object of the GenAnTK class:

The class is initialized with the paths to the following input parameters:
1. gatk_path (str) - Path to the GATK executable
2. out_path (str) - Path to the output directory
3. java (str) - Path to the java executable


Create an object:

```{python}
>>> from pyamd.gatk import GenAnTK
>>> gatk = GenAnTK('gatk', 'tests/variantcalling', 'java')
```

## Calling SNPs using HaplotypeCaller from GATK

Given a BAM file, reference fasta and sample name, the hapCaller module can be
used to call SNPs. The module return the path to the VCF file, and the
returncode for the process execution. The wrapper takes three input parameters:
1. bam_path (str) - Path to the read group annotated and index BAM file
2. ref_path (str) - Path to the reference fasta file with index and sequence dictionary
3. sam_name (str) - Name of the sample

```{python}
>>> vcf_path, returncode = gatk.hapCaller('tests/variantcalling/output_FM_SR_DD_RG.bam', 'tests/variantcalling/ref.fa', 'output')
```

# Variant calling  using Samtools

Listed below are the test cases for the different modules in the Samtools class

## Creating an object of the Samtools class:

The class is initialized with the paths to the following input parameters:
1. sam_path (str) - Path to the Samtools executable
2. bft_path (str) - Path to the Bcftools executable
3. out_path (str) - Path to the output directory


Create an object:

```{python}
>>> from pyamd.samtools import Samtools
>>> samtools = Samtools('samtools', 'bcftools', 'tests/variantcalling')
```

## Creating a pileup from BAM files:

Given a BAM file, reference fasta and sample name, the pileup module can be
used to create a mpileup for the sample. The module returns the path to a bcf
file and the returncode for the execution. the input parameters for the wrapper
are:
1. ref_path (str) - Path to the indexed fasta file
2. bam_path (str) - Path to the read group annotated BAM file
3. sam_name (str) - Sample name

```{python}
>>> bcf_path, returncode = samtools.pileup('tests/variantcalling/ref.fa', 'tests/variantcalling/output_FM_SR_DD_RG.bam', 'output')
>>> samtools.bcfindex(bcf_path)
```

## Creating variant calls from pileup BCF file

Given a BCF file, the bcftools module will generate SNP calls for the sample.
The module returns a VCF file and the returncode for the execution. The input
parameters for the wrapper are :
1. bcf_path (str) - Path to the indexed pileup BCF file

```{python}
>>> vcf_path, returncode = samtools.bcftools(bcf_path, 'output')
```

## Annotating VCF files with gene boundry information

NeST comes packaged with a custom VCF annotation tool that can annotate the Vcf
file with gene boundary, codon, amino acid change infomration. To create an
object of the annotater

```{python}
>>> from pyamd.parsers.vcf import Vcf
>>> annotater = Vcf.Annotate()
```

To annotate a VCF the following input parameters are required:
1. bed_path (str) - Path to the 12 column bed file for the reference genome
2. vcf_path (str) - Path to the input VCF file
3. fasta_path (str) - Path to the reference fasta file
4. out_path (str) - Path to the output directory

```{python}
>>> vcf_path = annotater.getAnnotation('tests/variantcalling/ref.bed', 'tests/variantcalling/output_variants_gatk.vcf', 'tests/variantcalling/ref.fa', 'tests/variantcalling')
```

## Merging variant calls from two different variant callers

NeST comes packaged with a custom VCF merging tool that can be used to merge any
two Vcf files for a given sample and create a consensus variant call. To create
an object for the merge tool the following input parameters are required:
1. fone (obj) - VCF reader object for the first VCF file
2. ftwo (obj) - VCF reader object for the second VCF file  
3. out_path(str) - Path to the output directory

```{python}
>>> from pyamd.parsers.vcf import Vcf
>>> fone = Vcf.Reader('tests/variantcalling/output_variants_gatk.vcf')
>>> ftwo = Vcf.Reader('tests/variantcalling/output_variants_samtools.vcf')
>>> merger = Vcf.Merge(fone, ftwo, 'tests/variantcalling').merge()
```
