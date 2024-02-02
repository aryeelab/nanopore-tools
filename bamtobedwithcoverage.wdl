version 1.0

workflow bamtobigwig {

    input {
        File genome
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
    call tobedgraph {
        input:
            genome = genome,
            filteredbai = filter.filteredbai,
            filteredbam = filter.filteredbam
    }
    output {
        File FivemCcpgbedgraph = tobedgraph.FivemCcpgbedgraph
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
    samtools import -T "*" ~{fastq} > temp.bam
    samtools fastq -T "*" temp.bam | minimap2 -ax map-ont -y ~{genome} --MD - | samtools sort -T tmp -o sorted.bam
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
task tobedgraph {
    input {
        File genome
        File filteredbam
        File filteredbai
    }
    command <<<
    modkit pileup ~{filteredbam} big5mC.cpg.bed --cpg --ref ~{genome} --ignore h -t 12 --combine-strands
    awk '$10 > 0 {printf "%s\t%d\t%d\t%d\t%2.3f\n" , $1,$2,$3,$10,$11}' big5mC.cpg.bed | sort -k1,1 -k2,2n > 5mC.cpg.bedgraph
    >>>
    runtime {
        docker: "us-central1-docker.pkg.dev/aryeelab/docker/modkit:latest"
		memory: "64G"
		disks: "local-disk 500 SSD"
		cpu: 8
    }
    output {
        File FivemCcpgbedgraph = "5mC.cpg.bedgraph"
    }
}
