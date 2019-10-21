library(readr, lib.loc="/usr/local/lib/R/site-library")
library(dplyr, lib.loc="/usr/local/lib/R/site-library")
library(stringr, lib.loc="/usr/local/lib/R/site-library")
library(GenomicRanges,lib.loc="/usr/local/lib/R/site-library")

nanopolish_out_kmer_in <- commandArgs(TRUE)[1]
outfile_by_read <- commandArgs(TRUE)[2]
outfile_by_kmer <- commandArgs(TRUE)[3]
threshold_log_likelihood <- commandArgs(TRUE)[4] #2.5
cpg_islands_file <- commandArgs(TRUE)[5] #hg38_CpG_islands.bed
valid_chrs_file <- commandArgs(TRUE)[6] #"hg38_chrs_noM.txt"
b_compartment_file <- commandArgs(TRUE)[7] #mergedBlocks_hg38_Hansen_2011.bed"

labelname <- unlist(lapply(strsplit(basename(nanopolish_out_kmer_in), split="__"), "[[", 1))


cpg_islands <- read_tsv(cpg_islands_file, col_names = c("chr", "start", "end", "cpg_island_name"))
gr_cpg_islands <- GRanges(cpg_islands$chr, IRanges(cpg_islands$start, cpg_islands$end), "cpg_island_name" = cpg_islands$cpg_island_name)

valid_chrs <- read_lines(valid_chrs_file)

b_compartment <- read_tsv(b_compartment_file, col_names = c("chr", "start", "end", "bin", "other","other"))
b_compartment <- b_compartment %>% filter(chr %in% valid_chrs)

gr_b_compartment <- GRanges(b_compartment$chr, IRanges(b_compartment$start, b_compartment$end), "bin" = b_compartment$bin)

kmer_input <- read_tsv(nanopolish_out_kmer_in,  col_types = "cciicdddiic", guess_max = 2147483647)

kmer_input <- kmer_input %>% filter(paste0("chr",chromosome) %in% valid_chrs)

kmer_input$num_motifs_WCGW <- str_count(kmer_input$sequence, pattern = "[A,T]CG[A,T]")

gr_kmer_input <- GRanges(paste0("chr",kmer_input$chromosome), IRanges(kmer_input$start, kmer_input$end))

kmer_input$intersects_cpg_island <- countOverlaps(gr_kmer_input, gr_cpg_islands) > 0
kmer_input$intersects_b_compartment <- countOverlaps(gr_kmer_input, gr_b_compartment) > 0
kmer_input$all_WCGWs <- kmer_input$num_motifs == kmer_input$num_motifs_WCGW
kmer_input$label <- labelname

write_tsv(kmer_input, path=outfile_by_kmer)

kmer_input %>%
  group_by(read_name, chromosome, label) %>%
  summarize(first_motif_occ = min(start),
            last_motif_occ = max(end),
            read_length = max(end) - min(start),
            total_CGs = sum(num_motifs),
            total_WCGWs = sum(num_motifs_WCGW),
            total_kmers = n(),
            fraction_kmers_intersect_CpG_island = sum(intersects_cpg_island)/n(),
            fraction_kmers_intersect_b_compartment = sum(intersects_b_compartment)/n(),
            mean_log_likelihood = mean(log_lik_ratio),
            #fraction_kmers_hypomethylated = sum(log_lik_ratio < -1*threshold_log_likelihood)/sum(abs(log_lik_ratio) > threshold_log_likelihood),
            fraction_kmers_methylated = sum(log_lik_ratio > threshold_log_likelihood)/sum(abs(log_lik_ratio) > threshold_log_likelihood),
            #get fraction methylation in B compartment, fraction methylation in non-B compartment
            fraction_kmers_methylated_in_B_compkmers  = sum(log_lik_ratio[intersects_b_compartment] > threshold_log_likelihood)/sum(abs(log_lik_ratio)[intersects_b_compartment] > threshold_log_likelihood),
            fraction_kmers_methylated_in_nonB_compkmers  = sum(log_lik_ratio[!intersects_b_compartment] > threshold_log_likelihood)/sum(abs(log_lik_ratio)[!intersects_b_compartment] > threshold_log_likelihood)
  ) -> summary_by_read
summary_by_read$is_vlong_read <- summary_by_read$read_length > 1e6 
summary_by_read$is_high_WCGW <- summary_by_read$total_WCGWs  > 50


write_tsv(summary_by_read, path=outfile_by_read)

write_tsv(summary_by_read %>% ungroup() %>% select(chromosome, first_motif_occ, last_motif_occ, fraction_kmers_methylated), 
          path=paste0(outfile_by_read,".reads.bedGraph"), col_names = FALSE)
