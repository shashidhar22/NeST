# Next-generation Sequence-analysis Toolkit (NeST) : A standardized bioinformatics framework for analyzing SNPs in next-generation sequencing data

Advancements in next-generation sequencing have led to the development of numerous bioinformatics tools and pipelines. Current tools for variant calling offer high-quality solutions; however, many tools are tailored for model organisms. Here, we present NeST, a consensus-based variant calling tool with a plug-and-play framework for use with any organism with minimal user input. NeST consists of four modules, integrating open-source bioinformatics tools and a custom VCF parser, to generate high-quality consensus variant calls. NeST was validated using targeted-amplicon deep sequencing data from 245 Plasmodium falciparum isolates to identify single-nucleotide polymorphisms conferring drug resistance. NeST offers a light-weight pipeline for variant calling with standardized outputs and minimal computational demands for easy deployment for use with various organisms and applications. The following document outlines details of installation, results from MaRS dataset and usage of individual modules for analysis.

1. [Overview of NeST framework](#Overview)
2. [Availability of code and installation](#Installation)
3. [Your first analysis](#First)
4. [Input standardization](#inputs)
5. [NeST class structure](#classes)

<a id="Overview"></a>
## Overview of the NeST framework:

NeST is a python based modular framework for consensus based variant calling. The overall analysis framework is broken down into four major blocks.
1. PrepInputs
2. VarCallEngine
3. VCFToolkit
4. Summarize

![NeST framework overview](images/Kookaburra.png)

The figure outlines the four key blocks of NeST and the steps performed by each step. VarCallEngine and VCFToolkit are spawned in parallel for each sample that is being analyzed in the study. By default, 4 parallel threads are spawned, to account for minimum available computational resource. This can be altered as per availability of resources.

<a id="Installation"></a>
## Availability of code and installation:

1. Download git repository:

   Clone the master branch of this repository.
   ```
   git clone https://github.com/shashidhar22/NeST
   ```

2. Installation:

   NeST comes with a install script that can be run to setup miniconda and create the virtual environment required to run NeST. To setup up miniconda and the NeST virtual environment, run the following command from the NeST directory.

   ```
   ./install.sh
   ```

   If you already have Anaconda or Miniconda installed on your system. You can skip the install step and just create a new environment for NeST using the configuration file provided in ```lib/``` directory. To create the NeST environment, run the following command from the NeST directory.

   ```
   conda env create -n nest -f lib/nest_env.yaml
   ```

   Once the environment is create it can be activated using the command.

   ```
   conda activate nest
   ```

   To deactivate the environment just type.

   ```
   conda deactivate
   ```

   Note: NeST virtual environment currently uses `bioconda` and `conda-forge` channels. There are known conflicts with the conda `defaults` channel and `conda-forge` channel, which can lead to errors in creation of the environment. To overcome the issue, the NeST environment ignores the `defaults` channel. If you plan to modify the NeST environment, please update the configuration file provided or make sure to ignore the `defaults` channel.

<a id="First"></a>
## Your first analysis

   NeST was conceptualized to identify mutations that confer anti-malarial drug resistance in *P.falciparum* (Talundzic et al., 2018). It was also applied for the detection of antibiotic drug resistance in *M.tubercolosis* (Colman et al., 2015). To make it easier to recreate these studies, we have included script to execute NeST on these two datasets.

   1. MaRS on NeST:

      To re-create the MaRS study (Talundzic et al., 2018), run the following command from the NeST directory.

      ```
      ./runPF.sh
      ```

      This downloads 10 sample fastq files from the [MaRS bioproject]((https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA428490) and executes NeST on the samples. The results are stored under `local/MaRS_test` folder within the NeST directory. The 10 target amplicon sequencing datasets, serve as a perfect test case for NeST and the analysis should take about 10 minutes on any laptop*.

      *\*Depending on internet speed*

   2. Detecting mutations conferring drug resistance in *M.tuberculosis* clinical samples from Colman et al., 2015:

      To re-create the analysis from Colman et al., 2015, run the following command from the NeST directory.

      ```
      ./runTB.sh
      ```

      This will download 57 paired fastq files from the bioproject [PRJNA271805](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA271805), merges the different runs and executes NeST on the clinical samples from the paper. The results are stored under `local/ColmanEtAl` folder within the NeST directory. This analysis takes significantly longer than the MaRS study and requires around 80GB of disk space, due to the large amount of data that is download. Make sure you have a good internet connection and adequate amount of coffee before starting this study.

   3. Executing your own analysis using NeST:

      NeST can be executed on your own dataset using the following command:

      ```
      python3 nest.py -i <path to input directory/ sra accession list> -a <adapter fasta file> -r < reference fasta file> -b < reference bed file> -o <output directory path> --varofint <CSV file with variants of interest>
      ```

      NeST can be run without a variant of interest file, but it will not produce any summary figures. To get a list of options that can be used with NeST just type

      ```
      python3 nest.py -h
      ```

      The details about the required input formats are listed in the next section.

<a id="inputs"></a>
## Input standardization:

NeST is designed to reduce the amount of user intervention with regards to inputs that the user needs to provide. However to enable standardization of inputs across all organisms we require that a particular file format be followed for the three inputs listed below:

1. Fastq files:

   The PrepInputs module in NeST highly simplifies the management of fastq files. The module accepts two input formats.
   1. Input directory path:

      This just requires the user to provide the path to a folder containing fastq files. The files are recognized by the file extension, so the files must have either ```fq```, ```fq.gz```, ```fastq``` or ```fastq.gz``` file extensions. The name convention of paired file can be ```_1```, ```_r1```, or ```_R1```.

   2. SRA accession list:

      This list requires a ```.txt``` with a list of SRA experiments, with one SRA number per line. This can be export from the SRA run selector tool.
      An example SRA accession is provided under ```fq/MaRS_test/SRA_Acc_list.txt```.

2. BED format:

   The BED (Browser Extensible Data) is an easy and lightweight format to list annotations for a genome. NeST uses a full BED or BED 12 column format file as a guide to annotate variants with codon and amino acid changes. The example file listed below shows the details of how to structure the BED file. The separation of contig, gene and exon level information makes this format highly portable across genomes. The BED 12 column format for most organisms can be export from the [UCSC table browser](https://genome.ucsc.edu/cgi-bin/hgTables). A detail explanation of the BED format can be found [here](https://genome.ucsc.edu/FAQ/FAQformat.html#format1)

   ```
   #contig start stop gene score strand  CDSstart  CDSstop rbg NoOfExons ExonLengths ExonStarts
   PfCRT 1	3399	PfCRT	.	+	95	3191	0	13	90,268,172,132,71,75,82,50,56,92,44,54,76,	96,364,812,1157,1443,1638,1810,2020,2208,2413,2699,2891,3115,
   MT	1	5967	COXIII	.	-	734	1573	0	1	839,	734,
   MT	1	5967	COL	.	+	1933	3471	0	1	1538,	1933,
   MT	1	5967	CYTOb	.	+	3492	4622	0	1	1130,	3492,
   PfDHFR	1	1827	PfDHFR	.	+	1	1827	0	1	1827,	1,
   PfDHPS	1	2417	PfDHPS	.	+	1	2417	0	3	135,1868,115,	1,313,2302,
   PfK13	1	2181	PfK13	.	+	1	2181	0	1	2181,	1,
   PfMDR1	1	4260	PfMDR1	.	+	1	4260	0	1	4260,	1,
   ```

3. Variant of Interest:

   The Summarize module within NeST, allows for easy summarization of variants called from all samples in a study. If a user specifies a list of variants of interest, a separate table will be created for these set of variants. The variants can be specified in ```.tsv```, ```.csv```, ```.xlsx``` format. And follows the format listed below

   | Chrom  | Gene   | RefAA | AAPos | AltAA |
   |:------:|:------:|:-----:|:-----:|:-----:|
   | PfCRT  | PfCRT  |   C   |   72  |   S   |
   | PfCRT  | PfCRT  |   V   |   73  |   V   |
   | PfMDR1 | PfMDR1 |   N   |   86  |   Y   |
   | PfMDR1 | PfMDR1 |   Y   |   184 |   F   |
   | MT     | CYTOb  |   I   |   258 |   M   |
