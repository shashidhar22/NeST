#!/usr/bin/Rscript

#Load each library
library(ggplot2) # used to create plots
library(reshape2) # used for data transformation
library(tidyr) # used for data transformation
library(readr) # used to read in csv files into tibbles
library(dplyr) # used for data transformation
library(optparse) # used to designate flags
library(RColorBrewer)
# set flags
option_list <- list(
  make_option(c("-i", "--input_file"), type = "character", default = NULL,
              help = "file path to Study_depth.csv", metavar = "character"),
  make_option(c("-o", "--output_file"), type = "character", default = NULL,
              help = "file path to output figure image file", metavar = "character")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
if(is.null(opt$input_file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file for Study_depth.csv).n", call.=FALSE)
}
if(is.null(opt$output_file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (output file path).n", call.=FALSE)
}

#Input data from MaRS output
tbl <- read_csv(opt$input_file, col_names = T) # Read in Study_depth.csv data into a tibble (a type of data frame)

#Reformat data to long from wide format to inorder to make plots
tbl.long <- melt(tbl, id.vars = "Variant", variable.name = "Sample", value.name = "Readdepth")
tbl.long <- separate(tbl.long, Variant, c("Loci", "SNP"), sep = ":", remove = F)
tbl.long$CodonPos <- as.numeric(substr(tbl.long$SNP, 2, nchar(tbl.long$SNP)-1))
tbl.long$Variant <- factor(tbl.long$Variant, unique(tbl.long$Variant))
tbl.long <- arrange(tbl.long, Loci, CodonPos)

#Make boxplot
#Set desired colors useing hexidecimal values for your plot (6 colors for 6 genes/targets)
#cbPalette <- c("#FF0000", "#0000FF", "#009900", "#FFFF00", "#9933FF", "#FF9900")
#Adding three colors for 9 gene/targets
if(length(unique(tbl.long$Loci)) <= 12){
	cbPalette <- brewer.pal(length(unique(tbl.long$Loci)), "Paired")
} else {
	cbPalette <- brewer.pal(length(unique(tbl.long$Loci)), "Spectral")
}
#cbPalette <- c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf','#999999')

#Generate boxplot
ggplot(tbl.long, aes(x=Variant, y=Readdepth, fill=Loci)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90)) + ylab("Read Depth") + theme(text = element_text(size=10), axis.text.x = element_text(size = 6)) + scale_fill_manual(values=cbPalette)

#Save plot as image file
ggsave(opt$output_file, plot = last_plot(), dpi = 300, width = 60, height = 30, units = c("cm"), limitsize = T)
