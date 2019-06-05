#Required packages: library(readr), library(stringr), library(optparse)

library(optparse) # used to designate flags

# set flags
option_list <- list(
  make_option(c("-i", "--input_file"), type = "character", default = NULL,
              help = "file path to Study_depth.csv", metavar = "character"),
  make_option(c("-r", "--report_snp"), type = "character", default = NULL,
              help = "file path to reportable snps list", metavar = "character"),
  make_option(c("-o", "--output_dir"), type = "character", default = NULL,
              help = "file path to output directory", metavar = "character")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
if(is.null(opt$input_file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file for Study_depth.csv)", call.=FALSE)
}
if(is.null(opt$report_snp)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (file path to reportable_SNPs.csv)", call. = FALSE)
}
if(is.null(opt$output_dir)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (output dir path)", call.=FALSE)
}
# Load readr library
library(readr)
library(stringr)


#Move reading of reportable SNPs to before changing the working directory
reportable_SNPs = read_csv(opt$report_snp, col_names = T)


###### read in summarized vcf data

vcfdata <- read_csv(opt$input_file, col_names = T)

## Set working directory
setwd(opt$output_dir)

## recode DP as numeric variable
vcfdata$DP[vcfdata$DP == "-inf"] = "-Inf"
vcfdata$DP = as.numeric(vcfdata$DP)

### NA and depth > 1 (logDP >= 0) means no mutation found
vcfdata$AF[is.na(vcfdata$AF) & vcfdata$DP >= 0] = 0
vcfdata$AF[vcfdata$AF == "" & vcfdata$DP >= 0] = 0


################ minimum DP?
#vcfdata$AF[vcfdata$DP == 0] = NA
vcfdata$AF[vcfdata$DP < 5] = NA


### convert frequency to proportion
vcfdata$AF = as.numeric(vcfdata$AF)/100

sampleIDs = unique(vcfdata$Sample)


###### read in list of reportable SNPs
#reportable_SNPs = read_csv(opt$report_snp, col_names = T)
reportable_SNPs = reportable_SNPs[reportable_SNPs$Gene != "",]
reportable_SNPs = cbind(reportable_SNPs,fullname = paste(reportable_SNPs$Gene,":", reportable_SNPs$RefAA, reportable_SNPs$AAPos, reportable_SNPs$AltAA ,sep=""))
genes = unique(reportable_SNPs$Gene)

head(reportable_SNPs)

process_vcf = function(sampleid) {
  ###### import

  vcftable = subset(vcfdata, Sample == sampleid)
  frequencies = sapply(1:dim(reportable_SNPs)[1], function (x) vcftable$AF[as.character(vcftable$Variant) == as.character(reportable_SNPs$fullname[x])])
  frequencies_vector = rep(NA, length(frequencies))
  frequencies_vector[unlist(lapply(frequencies,length)) == 1] = unlist(frequencies[unlist(lapply(frequencies,length)) == 1] )

  other_mutations = as.data.frame(cbind(fullname = as.character(vcftable$Variant)[!(as.character(vcftable$fullname) %in% as.character(reportable_SNPs$fullname))],
                                        AlFreq = vcftable$AF[!(as.character(vcftable$Variant) %in% as.character(reportable_SNPs$fullname))]))

  list(frequencies_vector,other_mutations)
}


processed_vcfs = sapply(sampleIDs,process_vcf)
frequencies_matrix = do.call(cbind,processed_vcfs[1,])

##### Table 1: Report Presence/Absence of Known Antimalarial Resistance SNPs

table1 = matrix("",length(reportable_SNPs$fullname),5)
table1[,1] = as.character(reportable_SNPs$Gene)
table1[,2] = as.character(paste(reportable_SNPs$RefAA,reportable_SNPs$AAPos, reportable_SNPs$AltAA,sep=""))
table1[,3] = paste(rowSums(frequencies_matrix == 0,na.rm=TRUE), "/", rowSums(!is.na(frequencies_matrix == 0)), " (",format(100*rowSums(frequencies_matrix == 0,na.rm=TRUE)/rowSums(!is.na(frequencies_matrix == 0)),digits=0,trim=TRUE),")",sep="")
table1[,4] = paste(rowSums(frequencies_matrix >0 & frequencies_matrix < 0.5,na.rm=TRUE), "/", rowSums(!is.na(frequencies_matrix == 0)), " (",format(100*rowSums(frequencies_matrix >0 & frequencies_matrix < 0.5,na.rm=TRUE)/rowSums(!is.na(frequencies_matrix == 0)),digits=0,trim=TRUE),")",sep="")
table1[,5] = paste(rowSums(frequencies_matrix <=1 & frequencies_matrix >= 0.5,na.rm=TRUE), "/", rowSums(!is.na(frequencies_matrix == 0)), " (",format(100*rowSums(frequencies_matrix <=1 & frequencies_matrix >= 0.5,na.rm=TRUE)/rowSums(!is.na(frequencies_matrix == 0)),digits=0,trim=TRUE),")",sep="")
colnames(table1) = c("Gene","SNP", "All wild type", "Minor mutant allele", "Major mutant allele")
write.csv(table1, "Reportable_SNPs_Report.csv")

### bar graph of Table 1
pdf("Reportable_SNPs.pdf",width = 20,height = 30)
frequency_summary = cbind(rowSums(frequencies_matrix == 0,na.rm=TRUE),rowSums(frequencies_matrix >0 & frequencies_matrix < 0.5,na.rm=TRUE),rowSums(frequencies_matrix <=1 & frequencies_matrix >= 0.5,na.rm=TRUE))
frequency_summary = frequency_summary[seq(from = dim(frequency_summary)[1], to = 1, by = -1),]
temp = barplot(t(frequency_summary/rowSums(frequency_summary)),col=c("cyan","green","pink"),horiz=TRUE,las=1,xlab="Frequency",xpd=TRUE)
axis(2,at = temp,labels = rev(reportable_SNPs$fullname),las=2,cex.lab = 0.25,tick = FALSE,cex.axis=0.6,line = -1)
axis(4,at = temp,labels = paste("N=",rowSums(frequency_summary),sep=""),las=2,cex.lab = 0.25,tick = FALSE,cex.axis=0.6,line = -0.5)
legend(0.4,max(temp)*1.1,legend=c("All wild type", "Minor mutant allele", "Major mutant allele"),cex=0.75,bty="n",fill =c("cyan","green","pink"),ncol=3,xpd=TRUE)
dev.off()
