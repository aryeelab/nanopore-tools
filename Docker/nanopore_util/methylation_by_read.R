library(readr, lib.loc="/usr/local/lib/R/site-library")
library(dplyr, lib.loc="/usr/local/lib/R/site-library")
library(stringr, lib.loc="/usr/local/lib/R/site-library")


infile <-commandArgs(TRUE)[1] #"/data/aryee/sowmya/ctc_nanopore/test-run-1__barcode06.dedup.methylation_calls.tsv"
outfile_by_read <- commandArgs(TRUE)[2]
input <- read_tsv(infile)


input$log_lik_ratio_repeated <- as.character(lapply(mapply(rep, input$log_lik_ratio, input$num_motifs), paste, collapse=","))

input$num_methylated_2.5_thresh <- as.numeric(lapply(strsplit(input$log_lik_ratio_repeated, split = ","), get_count <- function(x, threshold) { return (sum(as.numeric(x) > threshold)) }, threshold = 2.5))
input$num_UNmethylated_2.5_thresh <- as.numeric(lapply(strsplit(input$log_lik_ratio_repeated, split = ","), get_count <- function(x, threshold) { return (sum(as.numeric(x) < threshold)) }, threshold = -2.5))

input$num_methylated_0_thresh <- as.numeric(lapply(strsplit(input$log_lik_ratio_repeated, split = ","), get_count <- function(x, threshold) { return (sum(as.numeric(x) > threshold)) }, threshold = 0))
input$num_UNmethylated_0_thresh <- as.numeric(lapply(strsplit(input$log_lik_ratio_repeated, split = ","), get_count <- function(x, threshold) { return (sum(as.numeric(x) < threshold)) }, threshold = 0))


input$num_motifs_WCGW <- str_count(input$sequence, pattern = "[A,T]CG[A,T]")
input$log_lik_ratio_repeated_WCGW <- as.character(lapply(mapply(rep, input$log_lik_ratio, input$num_motifs_WCGW), paste, collapse=","))

input$num_methylated_2.5_thresh_WCGW <- as.numeric(lapply(strsplit(input$log_lik_ratio_repeated_WCGW, split = ","), get_count <- function(x, threshold) { return (sum(as.numeric(x) > threshold)) }, threshold = 2.5))
input$num_UNmethylated_2.5_thresh_WCGW <- as.numeric(lapply(strsplit(input$log_lik_ratio_repeated_WCGW, split = ","), get_count <- function(x, threshold) { return (sum(as.numeric(x) < threshold)) }, threshold = -2.5))

input$num_methylated_0_thresh_WCGW <- as.numeric(lapply(strsplit(input$log_lik_ratio_repeated_WCGW, split = ","), get_count <- function(x, threshold) { return (sum(as.numeric(x) > threshold)) }, threshold = 0))
input$num_UNmethylated_0_thresh_WCGW <- as.numeric(lapply(strsplit(input$log_lik_ratio_repeated_WCGW, split = ","), get_count <- function(x, threshold) { return (sum(as.numeric(x) < threshold)) }, threshold = 0))



input %>%
  group_by(read_name, chromosome) %>%
  summarize(start = min(start),
            end = max(end),
#            log_lik_ratio_collapsed = paste(log_lik_ratio_repeated, collapse=","),
            total_CGs = sum(num_motifs),
            total_CGs_methylated_2.5_thresh = sum(num_methylated_2.5_thresh),
            total_CGs_methylated_0_thresh = sum(num_methylated_0_thresh),
	          fraction_CGs_methylated_2.5_thresh = sum(num_methylated_2.5_thresh)/sum(num_motifs),
      	    fraction_CGs_methylated_0_thresh = sum(num_methylated_0_thresh)/sum(num_motifs),
            total_WCGW = sum(num_motifs_WCGW),
            total_WCGW_methylated_2.5_thresh = sum(num_methylated_2.5_thresh_WCGW),
            total_WCGW_methylated_0_thresh = sum(num_methylated_0_thresh_WCGW),
            fraction_WCGW_methylated_2.5_thresh = sum(num_methylated_2.5_thresh_WCGW)/sum(num_motifs_WCGW),
            fraction_WCGW_methylated_0_thresh = sum(num_methylated_0_thresh_WCGW)/sum(num_motifs_WCGW)) %>%
  write_tsv(path=outfile_by_read)
