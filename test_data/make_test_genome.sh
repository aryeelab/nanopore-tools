bam1="../cromwell-executions/preprocess_flowcell/34e08956-f94d-4941-9eb9-e1770a02c216/call-align/shard-0/execution/test-run-1__barcode06.dedup.bam"
bam2="../cromwell-executions/preprocess_flowcell/34e08956-f94d-4941-9eb9-e1770a02c216/call-align/shard-1/execution/test-run-1__unclassified.dedup.bam"

samtools view -F4 $bam1 | awk '{OFS="\t";print $3, $4-10000, $4+90000}' > regions.bed
samtools view -F4 $bam2 | awk '{OFS="\t";print $3, $4-10000, $4+90000}' >> regions.bed

bedtools getfasta -fi Homo_sapiens.GRCh38.dna.primary_assembly.fa -bed regions.bed > test_genome.fa
