# R Code to generate figures using data outputs from NeST pipeline  

### Reportable SNPs
#### Visualize the distribution of read depths for each reportable variant at each NeST loci
1. Run R script DepthPerReportSNP.R
   * Rscript DepthPerReportSNP.R -i /file/path/to/Study_depth.csv -o /file/path/to/output/image/file.png
#### Visualize the frequency of wild-type, major, and minor alleles for each variant (if > 1 variant present)
1. Download reportable_SNPs.csv file from this repo
2. Run R script reportableSNPsFreq.R
   * Rscript reportableSNPsFreq.R -i /file/path/to/Study_variants.csv -r /file/path/to/reportable_SNPs.csv -o /path/to/ouput/directory/

### Novel synonymous and non-synonymous exonic and intronic SNPs
#### Visualize the frequency of major and minor alleles for each novel variant (if > 1 variant present)
1. Run R script NovelExonicNonSynSNPs.R
   * Rscript NovelExonicNonSynSNPs.R -i /file/path/to/Study_novel_exonic_variants.csv -o /path/to/output/directory/
2. Run R script NovelExonicSynSNPs.R
   * Rscript NovelExonicSynSNPs.R -i /file/path/to/Study_novel_exonic_variants.csv -o /path/to/output/directory/
3. Run R script NovelIntronicSNPs.R
   * Rscript NovelIntronicSNPs.R -i /file/path/to/Study_novel_intronic_variants.csv -o /path/to/output/directory/

### Contact Information

* Contact Sarah Schmedes (obf2@cdc.gov) with any questions.
