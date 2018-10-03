# Summarize test cases

Listed below are the test cases for the different modules in the summarize class

## Creating an object of the summary class:

The class is initialized with the paths to the following input files:
1. fasta (str) - Path to reference Fasta file. The Fasta file for the MaRS
project can be found in ```ref/mdr.fa```
2. bed (str) - Path to reference BED file. The order of the records in the BED
file should correspond to the order of contigs in the Fasta file, and the format
of the BED record should follow the 12 column BED format as specified in the
README.  The BED file for the MaRS project can be found in ```ref/mdr.bed```
3. voi (str) - Path to variant of interest file. The variant of interest file
can be in CSV, TSV, and XLSX format, with five columns, CHROM, GENE, RefAA,
AltAA and AAPos. The variant of interest file for the MaRS project can be found
in ```ref/Reportable_SNPs.csv```
4. out_path (str) - Path to study output directory. The directory structure for
the out_path folder is as follows
```
.
+-- kookaburra.log
+-- sample_directory
    +-- alignments
        +-- output.sam
    +-- CleanedFastq
        +-- sample_1.fastq_cleaned.fq
        +-- sample_2.fastq_cleaned.fq
        +-- sample_2.fastq_stats.txt
    +-- completion
        +-- align.rt
        +-- bbduk.rt
        +-- bcfcall.rt
        +-- bcfindex.rt
        +-- fixmate.rt
        +-- gatk.rt
        +-- pileup.rt
        +-- readgroup.rt
        +-- sort.rt
    +-- mod_fasta.fa
    +-- out_fixmate.bam
    +-- out_fixmate_sorted.bam
    +-- out_fixmate_sorted_RG.bai
    +-- out_fixmate_sorted_RG.bam
    +-- sample_variants_annotated.vcf
    +-- sample_variants_gatk_annotated.vcf
    +-- sample_variants_gatk.vcf
    +-- sample_variants_gatk.vcf.idx
    +-- sample_variants_merged_annotated.vcf
    +-- sample_variants.vcf
    +-- variants.bcf
+-- sample2_directory
    +--
```

Create an object:
```{python}
>>> from pyamd.prepinputs import Prepper
>>> from pyamd.summarize import Summary
>>> config = Prepper('tests/prepinputs/SRR_Acc_List.txt', None).prepInputs()
>>> summarizer = Summary('ref/pfalciparum/mdr.fa', 'ref/pfalciparum/mdr.bed', 'ref/pfalciparum/Reportable_SNPs.csv', 'tests/summarize')
```

## Get the codon base range for a given amino acid position

The `getBaseRange` module can be used to get the genomic base range of the codon
for a given amino acid position. The module returns a range with a 0-inclusive
start and 0-exclusive stop.

```{python}
>>> summarizer.getBaseRange('PfCRT', 'PfCRT', 271)
range(1710,1713)
```

## Generate a pandas dataframe from the variant of interest file

The `getVarOfInt` module can be used to generate a pandas dataframe from the
variant of interest file. The module accepts variants of interest file in
XLSX, CSV, and TSV format

```{python}
>>> voi_table = summarizer.getVarOfInt()
>>> voi_table.head()
            Chrom   Gene RefAA  AAPos AltAA  GenomicStart  GenomicStop   SNP
Variant                                                                     
PfCRT:C72S  PfCRT  PfCRT     C     72     S           NaN          NaN  C72S
PfCRT:V73V  PfCRT  PfCRT     V     73     V           NaN          NaN  V73V
PfCRT:M74I  PfCRT  PfCRT     M     74     I           NaN          NaN  M74I
PfCRT:N75E  PfCRT  PfCRT     N     75     E           NaN          NaN  N75E
PfCRT:K76T  PfCRT  PfCRT     K     76     T           NaN          NaN  K76T
```

## Get per gene coverage for a sample dataset

The `getGeneStats` module can be used to generate gene-wise dictionary with the
per-base coverage over the gene segment.

```{python}
>>> summarizer.getGeneStats('tests/summarize/SRR6463549/output_fixmate_sorted_RG.bam')
```

## Get per gene per base coverage for all samples in a study

The `getExpCoverage` module can be used to generate sample-wise, gene-wise
dictionary of dictionary of list containing the per base coverage of each gene
in each sample of a given study

```{python}
>>> coverageDict = summarizer.getExpCoverage()
```

## Identify which gene of which sample passed a threshold for coverage

The `checkDepthPass` module can be used to generate a dictionary of dictionary
describing whether the gene passed depth threshold per gene per study.

```{python}
>>> passDict  = summarizer.checkDepthPass()
```

## Generate a dataframe containing all intronic variant calls in a study

Given a study output path, generate a pandas dataframe containing
all intronic variants. The module gets a list of variants from the
output folder. A dictionary of list is initialized with the fields that
will be recorded in the intron summary tables. The features of the
intron summary tables are listed below:
    1. Chrom: Chromosome where the variant is found.
    2. Gene: Gene boundary within which the variant is found.
    3. Pos: Variant position.
    4. Qual: Phred based quality score of variant call.
    5. Ref: Reference base call for the variant call.
    6. Alt: Alternate(variant) base call for the variant call.
    7. Exon: Exon boundary within which the variant is found. This
             value is set to "Intron" for the intron table.
    8. AAPos: Amino acid position for the variant call. This value is
              set to np.nan value for the intron table.
    9. RefCodon: Reference triplet codon for the variant call. This
                 value is set to "NA" for the intron table.
   10. AltCodon: Alternate(variant) base call for the variant call. This
                 value is set to "NA" for the intron table.
   11. RefAA: Reference amino acid for the variant call. This value is
              set to "NA" for the intron table.
   12. AltAA: Alternate(variant) base call for the variant call. This
              value is set to "NA" for the intron table.
   13. DP: Depth for the codon to which the variant call belong.
   14. AF: Allele frequency for the variant call.
   15. Conf: Number of variant callers that identified the variant call.
The table is indexed by the following fields:
    1. Sample: Sample name
    2. Variant: A string containing the variant information in
                <Chrom>:<Ref><Pos><Alt> for the intron table.

```{python}
>>> intron_table = summarizer.getIntronTables()
>>> intron_table.head()
                         AAPos          AF Alt AltAA AltCodon  Chrom  Conf  \
Sample     Variant                                                           
SRR6463549 PfCRT:G211T     NaN  100.000000   T    NA       NA  PfCRT     2   
           PfCRT:T1079A    NaN   75.000000   A    NA       NA  PfCRT     1   
           PfCRT:T1083A    NaN  100.000000   A    NA       NA  PfCRT     1   
           PfCRT:A2193T    NaN   50.000000   T    NA       NA  PfCRT     2   
           MT:G1692A       NaN   99.150142   A    NA       NA     MT     2   

                          DP    Exon   Gene   Pos        Qual Ref RefAA  \
Sample     Variant                                                        
SRR6463549 PfCRT:G211T    53  Intron  PfCRT   211   1376.7700   G    NA   
           PfCRT:T1079A   26  Intron  PfCRT  1079     39.8786   T    NA   
           PfCRT:T1083A   31  Intron  PfCRT  1083     74.0000   T    NA   
           PfCRT:A2193T  114  Intron  PfCRT  2193    507.7700   A    NA   
           MT:G1692A     353  Intron     MT  1692  11298.7700   G    NA   

                        RefCodon  
Sample     Variant                
SRR6463549 PfCRT:G211T        NA  
           PfCRT:T1079A       NA  
           PfCRT:T1083A       NA  
           PfCRT:A2193T       NA  
           MT:G1692A          NA  
```


## Generate a dataframe containing all exonic variant calls in a study

Given a study output path, generate a pandas dataframe containing
all known and novel exonic variants. The module gets a list of VCF files
from the output folder. A dictionary of list is initialzed with the
fields that will be recorded in the intron summary tables. If a variant
from the list of variants of interest is not found in the sample, an
empty record is added to the dataframe with np.nan or "NA" values for
all field expect the variant descriptor fields. The features
of the exonic summary tables are listed below:
  1. Chrom: Chromosome where the variant is found.
  2. Gene: Gene boundary within which the variant is found.
  3. Pos: Variant position.
  4. Qual: Phred based quality score of variant call.
  5. Ref: Reference base call for the variant call.
  6. Alt: Alternate(variant) base call for the variant call.
  7. Exon: Exon boundary within which the variant is found. This value is set to exon number for the known and novel exonic table.
  8. AAPos: Amino acid position for the variant call.
  9. RefCodon: Reference triplet codon for the variant call.
  10. AltCodon: Alternate(variant) base call for the variant call.
  11. RefAA: Reference amino acid for the variant call.
  12. AltAA: Alternate(variant) base call for the variant call.
  13. DP: Depth for the codon to which the variant call belong.
  14. AF: Allele frequency for the variant call.
  15. Conf: Number of variant callers that identified the variant call.
The table is indexed by the following fields:
  1. Sample: Sample name
  2. Variant: A string containing the variant information in Chrom>:<RefAA><AAPos><AltAA> for the intron table.

```{python}
>>> exon_table = summarizer.getVarTables()
>>> exon_table = summarizer.getVarTables()
                         AAPos         AF  Alt AltAA AltCodon   Chrom  Conf  \
Sample     Variant                                                            
SRR6463549 MT:C18C          18   98.44358    C     C      TGC      MT     2   
           PfDHFR:N51I      51  100.00000    T     I      ATT  PfDHFR     2   
           PfDHFR:C59R      59  100.00000    C     R      CGT  PfDHFR     2   
           PfDHFR:S108N    108   99.87163    A     N      AAC  PfDHFR     2   
           PfDHPS:A437G    437  100.00000  G,T     G      GGT  PfDHPS     2   

                          DP   Exon    Gene   Pos      Qual  Ref RefAA  \
Sample     Variant                                                       
SRR6463549 MT:C18C       257  exon1   CYTOb  3545   8581.77    T     C   
           PfDHFR:N51I   741  exon1  PfDHFR   152  30550.77    A     N   
           PfDHFR:C59R   762  exon1  PfDHFR   175  31171.77    T     C   
           PfDHFR:S108N  792  exon1  PfDHFR   323  24018.77    G     S   
           PfDHPS:A437G  633  exon2  PfDHPS  1486  28251.77  C,G     A   

                        RefCodon  
Sample     Variant                
SRR6463549 MT:C18C           TGT  
           PfDHFR:N51I       AAT  
           PfDHFR:C59R       TGT  
           PfDHFR:S108N      AGC  
           PfDHPS:A437G      GCG  
```

## Create a table showing the presence or absence of variants of interest across all samples in a study

Returns a table containing sample wise breakdown about the presence
or absence of variants of interest in the study. If there are samples
with no calls for any particular variants of interest, if the sample has
coverage at the variant location and no variant call, record is labeled
as a "WT" call in the FinalCall field. If there is a variant call at any
of the locations, the record is labeled as "SNP" call in the FinalCall
field. This module internally calls the following modules:
  1. `getVarTables()`
  2. `getVarOfInt()`

The module produces a pandas dataframe with the following fields:
  1. Sample : Sample name
  2. Variant : Variant of interest
  3. Chrom : Chromosome
  4. Gene : Gene name
  5. SNP: SNP of interest
  6. FinalCall : Whether the sample had the SNP or it was a wild type call (WT)
  7. Ref : Reference genomic base
  8. Alt : Alternate(Variant) genomic base
  9. Pos : Genomic position of variant
  10. Qual : Phred based quality score
  11. RefCodon : Reference codon at variant of interest
  12. RefAA : Reference amino acid at variant of interest
  13. AltCodon : Alternate codon at variant of interest
  14. AltAA : Alternate amino acid at variant of interest
  15. AAPos : Amino acid position of interest
  16. Exon : Exon number where the variant is found
  17. AF : Allele frequency of variant of interest
  18. DP : Depth of coverage of variant of interest
  19. Conf : Number of variant callers that called the variant of interest        

```{python}
>>> knownVar_table = summarizer.getRepSnps()
>>> knownVar_table.head()
                       Chrom   Gene   SNP FinalCall  Ref  Alt  Pos  Qual  \
Sample     Variant                                                         
SRR6463549 PfCRT:C72S  PfCRT  PfCRT  C72S        WT  NaN  NaN  NaN   NaN   
           PfCRT:V73V  PfCRT  PfCRT  V73V        WT  NaN  NaN  NaN   NaN   
           PfCRT:M74I  PfCRT  PfCRT  M74I        WT  NaN  NaN  NaN   NaN   
           PfCRT:N75E  PfCRT  PfCRT  N75E        WT  NaN  NaN  NaN   NaN   
           PfCRT:K76T  PfCRT  PfCRT  K76T        WT  NaN  NaN  NaN   NaN   

                      RefCodon RefAA AltCodon AltAA  AAPos Exon  AF  DP  
Sample     Variant                                                       
SRR6463549 PfCRT:C72S      NaN     C      NaN     S   72.0  NaN NaN NaN  
           PfCRT:V73V      NaN     V      NaN     V   73.0  NaN NaN NaN  
           PfCRT:M74I      NaN     M      NaN     I   74.0  NaN NaN NaN  
           PfCRT:N75E      NaN     N      NaN     E   75.0  NaN NaN NaN  
           PfCRT:K76T      NaN     K      NaN     T   76.0  NaN NaN NaN  
```

## Create a table of novel exonic and intronic variants across all samples in a study

Returns a table containing sample wise breakdown about the presence
of novel variants in the study. This table will contain a FinalCall
field, since only variants are reported and not wild type calls. This
module internally calls the following modules:
  1. `getVarTables()`
  2. `getVarOfInt()`

The module produces a pandas dataframe with the following fields:
  1. Sample : Sample name
  2. Variant : Variant of interest
  3. Chrom : Chromosome
  4. Gene : Gene name
  5. Ref : Reference genomic base
  6. Alt : Alternate(Variant) genomic base
  7. Pos : Genomic position of variant
  8. Qual : Phred based quality score
  9. RefCodon : Reference codon at variant of interest
  10. RefAA : Reference amino acid at variant of interest
  11. AltCodon : Alternate codon at variant of interest
  12. AltAA : Alternate amino acid at variant of interest
  13. AAPos : Amino acid position of interest
  14. Exon : Exon number where the variant is found
  15. AF : Allele frequency of variant of interest
  16. DP : Depth of coverage of variant of interest
  17. Conf : Number of variant callers that called the variant of interest

```{python}
>>> novelVar_table = summarizer.getNovSnps()
>>> novelVar_table.head()
                   Chrom   Gene Ref Alt   Pos     Qual RefCodon RefAA  \
Sample     Variant                                                      
SRR6463549 MT:C18C    MT  CYTOb   T   C  3545  8581.77      TGT     C   

                   AltCodon AltAA  AAPos   Exon        AF   DP  Conf  
Sample     Variant                                                    
SRR6463549 MT:C18C      TGC     C     18  exon1  98.44358  257     2  
```

## Get average coverage across a chromosome in sample

The `getBamStat` module can be used to get the average coverage across a chromosome
in a given sample.

```{python}
>>> bamStat = summarizer.getBamStat('tests/summarize/SRR6463549/output_fixmate_sorted_RG.bam', 'PfCRT', 1, 3399)
>>> bamStat
3524
```

## Annotate a variant with the average codon depth, rather than the base position

Given a study dataframe, annotate all the variants with the average triplet
depth, rather than the using just the depth at the base as a metric. The
resulting dataframe has the same fields as the `getRepSnps` dataframe

```{python}
>>> kvDepth_table = summarizer.getDepthStats(knownVar_table)
                       Chrom   Gene   SNP FinalCall  Ref  Alt  Pos  Qual  \
Sample     Variant                                                         
SRR6463549 PfCRT:C72S  PfCRT  PfCRT  C72S        WT  NaN  NaN  NaN   NaN   
           PfCRT:V73V  PfCRT  PfCRT  V73V        WT  NaN  NaN  NaN   NaN   
           PfCRT:M74I  PfCRT  PfCRT  M74I        WT  NaN  NaN  NaN   NaN   
           PfCRT:N75E  PfCRT  PfCRT  N75E        WT  NaN  NaN  NaN   NaN   
           PfCRT:K76T  PfCRT  PfCRT  K76T        WT  NaN  NaN  NaN   NaN   

                      RefCodon RefAA AltCodon AltAA  AAPos Exon  AF   DP  Conf  
Sample     Variant                                                              
SRR6463549 PfCRT:C72S      NaN     C      NaN     S   72.0  NaN NaN  135   2.0  
           PfCRT:V73V      NaN     V      NaN     V   73.0  NaN NaN  136   2.0  
           PfCRT:M74I      NaN     M      NaN     I   74.0  NaN NaN  135   2.0  
           PfCRT:N75E      NaN     N      NaN     E   75.0  NaN NaN  134   2.0  
           PfCRT:K76T      NaN     K      NaN     T   76.0  NaN NaN  134   2.0  
```

## Create a study JSON that summarizes the whole study in a structured, portable manner

Create a summary compressed JSON file for all variants in all samples
in the study. The information recorded in the JSON format include:
  1. Study name
  2. Number of samples in study
  3. Date of analysis
  4. Sample name:
    1. Variant description string in <Gene><RefAA><AAPos><AltAA> format:
      1. Reference amino acid
      2. Alternate amino acid
      3. Variant call position
      4. Final call field, indicating the wild type or variant call based of BAM file coverage
      5. Allele frequency for the variant call
      6. Depth for the variant call
      7. Confidence for the vairant call based on the number of variant callers that identified the variant


```{python}
>>> json_file = summarizer.toJSON()
```

## Create study CSV files separating the study variant calls in Known, Novel exonic and Novel intronic

Create CSV files from the DataFrames and generate the allele frequency and depth views for known and novel files. The `toCSV` module calls, `getRepSnps` module internally and generates the dataframe along with the CSV file. The module creates the following files :
  1. Study_known_variants.csv : Known variant calls
  2. Study_novel_exonic_variants.csv : Novel exonic variant calls
  3. Study_novel_intronic_variants.csv : Novel intronic variant calls
  4. Study_**var_type**\_variants_allele_frequency.csv : Allele frequency view for the study for all samples in the study, across all variants of that type
  5. Study_**var_type**\_variants_depth.csv : Depth of coverage view for the study for all samples in the study, across all variants of that type

```{python}
>>> summarizer.toCSV('known')
>>> summarizer.toCSV('novel')
```

## Summarize a study. Create tables, JSONs and summary figures for a study

The `getSummary` module is the main module for the Summary class, it internally calls the toCSV and toJSON modules to create the study summary tables and the study JSON. Apart from that, it creates the following figures:
  1. Study_depth.pdf : Summary figure showing depth of coverage for each of the variants of interest
  2. Reportable_SNPs.pdf : Summary figure showing the allele frequency of the variants of interest in the study sample set, showing the prevalence of the variant of interest in the study population
  3. Novel_SNPs_exonic_nonsyn.pdf : Summary figure showing the allele frequency of the novel exonic non-synonymous in the study sample set, showing the prevalence of these mutations in the study population
  4. Novel_SNPs_exonic_syn.pdf : Summary figure showing the allele frequency of the novel exonic synonymous in the study sample set, showing the prevalence of these mutations in the study population
  5. novel_SNPs_intronic.pdf : Summary figure showing the allele frequency of the novel intronic in the study sample set, showing the prevalence of these mutations in the study population

```{python}
>>> summarizer.getSummary()
```
