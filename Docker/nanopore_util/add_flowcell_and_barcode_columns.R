#!/usr/bin/env Rscript

library(stringr)

args <- commandArgs(trailingOnly=TRUE)
infile <- args[1]
outfile <- args[2]

if (is.na(infile)) 
  stop("Input file (argument 1) not specified")

if (is.na(outfile)) 
  stop("Output file (argument 2) not specified")


df <- read.csv(infile, stringsAsFactors = FALSE)
if ("sample_id" %in% colnames(df)) 
    stop("The input file already has a sample_id column")


if (nrow(df) > 0) {
    df <- df[,colnames(df) != "X"]
    x <- sub(".fq.gz", "", basename(df$fastq_gz))
    flowcell_and_barcode <- str_split_fixed(x, pattern = "__", n=2)
    colnames(flowcell_and_barcode) <- c("flowcell_id", "barcode_id")
    out <- cbind(sample_id="", flowcell_and_barcode, df)
    write.csv(out, outfile, quote=FALSE, row.names=FALSE)
    message("Wrote ", outfile)
} else {
    out <- data.frame(sample_id=character(), flowcell_id=character(), barcode_id=character(), df)
    write.csv(out, outfile, quote=FALSE, row.names=FALSE)
    message("Wrote ", outfile, " but there are no samples in it!")
}

