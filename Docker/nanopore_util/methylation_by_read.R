library(readr)
library(dplyr)
library(stringr)


infile <-commandArgs(TRUE)[1] #"/data/aryee/sowmya/ctc_nanopore/test-run-1__barcode06.dedup.methylation_calls.tsv"
outfile_by_read <- commandArgs(TRUE)[2]
input <- read_tsv(infile, col_types=
    cols(
      chromosome = col_character(),
      strand = col_character(),
      start = col_double(),
      end = col_double(),
      read_name = col_character(),
      log_lik_ratio = col_double(),
      log_lik_methylated = col_double(),
      log_lik_unmethylated = col_double(),
      num_calling_strands = col_double(),
      num_motifs = col_double(),
      sequence = col_character())
    )


if (nrow(input)==0) {
    message("Warning: No lines in the input file")
    file.create(outfile_by_read)
} else {
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

}

