version 1.0

workflow dorado_basecall {

    input {
        File genome
        File chromsizes
        String sample_id
        File fast5_archive
        String basecall_model
    }

    call basecall {
        input:
            sample_id = sample_id,
            fast5_archive = fast5_archive,
            basecall_model = basecall_model
    }
    call minimapalign {
        input:
            unmapped_bam = basecall.unmapped_bam,
            genome = genome
    }
    call filter {
        input:
            sortedbam = minimapalign.sortedbam
    }
    call bamtobigwig {
        input:
            chromsizes = chromsizes,
            genome = genome,
            filteredbai = filter.filteredbai,
            filteredbam = filter.filteredbam
    }
    output {
        File FivemCcpgbedgraph = bamtobigwig.FivemCcpgbedgraph
        File FivemCcpgbw = bamtobigwig.FivemCcpgbedgraph
        File unmapped_bam = basecall.unmapped_bam
        File? duplex_unmapped_bam = basecall.duplex_unmapped_bam
        File? pairs = basecall.pairs
    }

    meta {
        author: "Martin Aryee"
        email:"martin.aryee@gmail.com"
    }
}
task basecall  {
    input {
        String sample_id
        File fast5_archive
        String basecall_model
        Int disk_gb = ceil(size(fast5_archive, "GB")*3) + 5
    }
    command <<<
        filetype=$(file ~{fast5_archive})

        if [[ ~{fast5_archive} == *.pod5 ]]; then
            mkdir pod5s
            ln -s ~{fast5_archive} pod5s/reads.pod5
        fi

        if [[ "$filetype" == *"gzip compressed data"* ]]; then
          echo "FAST5s appear to be compressed with gzip. Decompressing..."
          mkdir fast5s
          tar zxvf ~{fast5_archive} -C fast5s
          # Convert FAST5 into POD5
          pod5 convert fast5 -r fast5s pod5s/reads.pod5 --threads 12 
        fi

        if [[ "$filetype" == *"Zip archive data"* ]]; then
          echo "FAST5s appear to be compressed with zip. Decompressing..."
          mkdir fast5s 
          unzip ~{fast5_archive} -d fast5s
          # Convert FAST5 into POD5
          pod5 convert fast5 -r fast5s pod5s/reads.pod5 --threads 12 
        fi
        
        # Simplex call with --emit-moves and --modified-bases
        dorado basecaller /dorado_models/~{basecall_model} pod5s --modified-bases 5mCG_5hmCG --emit-moves | samtools view -Sh > ~{sample_id}.unmapped.bam

        # Identify potential pairs
        duplex_tools pair --output_dir ./pairs ~{sample_id}.unmapped.bam
    
        # Stereo duplex basecall
        if [ -f "pairs/pair_ids_filtered.txt" ]; then
            NUM_PAIRS=$(wc -l pairs/pair_ids_filtered.txt)
            echo "${NUM_PAIRS} pairs found. Duplex calling..."
            dorado duplex /dorado_models/~{basecall_model} pod5s --pairs pairs/pair_ids_filtered.txt | samtools view -Sh > ~{sample_id}.duplex.unmapped.bam    
        else 
            echo "No pairs found."
        fi
    >>>
    runtime {
        gpuType: "nvidia-tesla-v100"
        gpuCount: 1
        cpu: 12
        disks: "local-disk " + disk_gb + " SSD" 
        memory: "32GB"
        nvidiaDriverVersion: "470.161.03"
        zones: ["us-central1-a"] 
        docker: "us-central1-docker.pkg.dev/aryeelab/docker/dorado"
    }
    output {
        File unmapped_bam = "~{sample_id}.unmapped.bam"
        File? duplex_unmapped_bam = "~{sample_id}.duplex.unmapped.bam"
        File? pairs = "pairs/pair_ids_filtered.txt"
    }
}
task minimapalign {
    input {
        File unmapped_bam
        File genome
    }
    command <<<
    samtools fastq -T "*" ~{unmapped_bam} | minimap2 -a -x map-ont ~{genome} | samtools sort -T tmp -o sorted.bam
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
task bamtobigwig {
    input {
        File genome
        File chromsizes
        File filteredbam
        File filteredbai
    }
    command <<<
    modbam2bed -t 12 -m 5mC --cpg ~{genome} ~{filteredbam} > big5mC.cpg.bed
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