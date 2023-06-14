version 1.0

workflow bamtobigwig {

    input {
        File genome
        File chromsizes
        File fastq
    }

    call minimapalign {
        input:
            fastq = fastq,
            genome = genome
    }
    call filter {
        input:
            sortedbam = minimapalign.sortedbam
    }
    call tobigwig {
        input:
            chromsizes = chromsizes,
            genome = genome,
            filteredbai = filter.filteredbai,
            filteredbam = filter.filteredbam
    }
    output {
        File FivemCcpgbedgraph = tobigwig.FivemCcpgbedgraph
        File FivemCcpgbw = tobigwig.FivemCcpgbw
    }

    meta {
        author: "Martin Aryee"
        email:"martin.aryee@gmail.com"
    }
}
task minimapalign {
    input {
        File fastq
        File genome
    }
    command <<<
    minimap2 -ax map-ont -y ~{genome} ~{fastq}| samtools sort -T tmp -o sorted.bam
    samtools index sorted.bam
    >>>
    runtime {
        docker: "us-central1-docker.pkg.dev/aryeelab/docker/minimap2:latest"
		memory: "64G"
		disks: "local-disk 500 SSD"
		cpu: 8
    }
    output {
        File sortedbam = "sorted.bam"
    }
}
task filter {
    input {
        File sortedbam
    }
    command <<<
    samtools view -bh -q 50 ~{sortedbam} > filtered.bam
    samtools index filtered.bam
    >>>
    runtime {
        docker: "us-central1-docker.pkg.dev/aryeelab/docker/samtools:latest"
		memory: "64G"
		disks: "local-disk 500 SSD"
		cpu: 8
    }
    output {
        File filteredbam = "filtered.bam"
        File filteredbai = "filtered.bam.bai"
    }
}
task tobigwig {
    input {
        File genome
        File chromsizes
        File filteredbam
        File filteredbai
    }
    command <<<
    modkit pileup ~{filteredbam} big5mC.cpg.bed --cpg --ref ~{genome} --ignore h -t 12 --combine-strands
    awk '$10 > 0 {printf "%s\t%d\t%d\t%2.3f\n" , $1,$2,$3,$11}' | sort -k1,1 -k2,2n > 5mC.cpg.bedgraph
    bedGraphToBigWig 5mC.cpg.bedgraph ~{chromsizes} 5mC.cpg.bw
    >>>
    runtime {
        docker: "us-central1-docker.pkg.dev/aryeelab/docker/bedtools:latest"
		memory: "64G"
		disks: "local-disk 500 SSD"
		cpu: 8
    }
    output {
        File FivemCcpgbw = "5mC.cpg.bw"
        File FivemCcpgbedgraph = "5mC.cpg.bedgraph"
    }
}