# Kookaburra : A standardized bioinformatics framework for analyzing SNPs in next-generation sequencing data

Advancements in next-generation sequencing have led to the development of numerous bioinformatics tools and pipelines. Current tools for variant calling offer high-quality solutions; however, many tools are tailored for model organisms. Here, we present Kookaburra, a consensus-based variant calling tool with a plug-and-play framework for use with any organism with minimal user input. Kookaburra consists of four modules, integrating open-source bioinformatics tools and a custom VCF parser, to generate high-quality consensus variant calls. Kookaburra was validated using targeted-amplicon deep sequencing data from 245 Plasmodium falciparum isolates to identify single-nucleotide polymorphisms conferring drug resistance. Kookaburra offers a light-weight pipeline for variant calling with standardized outputs and minimal computational demands for easy deployment for use with various organisms and applications. The following document outlines details of installation, results from MaRS dataset and usage of individual modules for analysis.

1. [Overview of Kookaburra framework](#Overview)
2. [Availability of code and installation](#Installation)
3. [Kookaburra for Malaria Resistance Surveillance(MaRS)](#MaRS)
4. [Input standardization](#inputs)
5. [Kookaburra class structure](#classes)

<a id="Overview"></a>
## Overview of the Kookaburra framework:

Kookaburra is a python based modular framework for consensus based variant calling. The overall analysis framework is broken down into four major blocks.
1. PrepInputs
2. VarCallEngine
3. VCFToolkit
4. Summarize

![Kookaburra framework overview](images/Kookaburra.png)

The figure outlines the four key blocks of Kookaburra and the steps performed by each step. VarCallEngine and VCFToolkit are spawned in parallel for each sample that is being analyzed in the study. By default, 4 parallel threads are spawned, to account for minimum available computational resource. This can be altered as per availability of resources.

<a id="Installation"></a>
## Availability of code and installation:

1. Download git repository:

Clone the master branch of this repository.
```{sh}
git clone https://github.com/shashidhar22/kookaburra
```

2. Setup virtualenv for Python3:

Kookaburra requires [Python3](https://www.python.org/downloads/) to be installed with [Pip](https://pip.pypa.io/en/stable/installing/) available. Please make sure this is available on the system.
To avoid clashes with system version of required python modules, we recommend using a virtualenv
Run the following command to install virtualenv, if you already have virtualenv installed

```{sh}
python3 -m pip install virtualenv
virtualenv kook_env                   # Setup mars virtual environment
source kook_env/bin/activate          # Activate virtual environment
```
> If successfully activated, you should see now (mars_venv) in front of your terminal username.

3. Installing python dependencies and downloading third party libraries:

Kookaburra uses many python modules
Run the following command to install the dependencies
```{sh}
pip3 install pysam matplotlib seaborn pandas numpy xlrd openpyxl
pip3 list --format=columns
```

4. Your first analysis:

Kookaburra comes packaged with an SRA accession list from the [MaRS]((https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA428490) experiment. The includes the SRA accession for 10 Illumina paired end samples. Running the command listed below, will download the 10 samples using [SRAToolkit](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software) and run kookaburra on it.
```{sh}
sh run.sh fq/MaRS_test/SRA_Acc.txt local/MaRs_test
```
To run Kookaburra on locally stored fastq files. You can just provide the path to the input directory instead of the accession list.
For example if you have stored your fastq files in ```fq/``` folder and you want to store the results in the folder ```local/```. You can run the following command from the Kookaburra directory.

```{sh}
sh run.sh fq/ local/
```
