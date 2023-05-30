version 1.0

workflow guppytonanopolish {

	input {
        File reads
        String type
        String typeagain
        File genome
        File config
        File model
        File sortedbed
        File pythonscript
        File chromsizes
	}

	call guppy {
		input:
			reads=reads,
            model=model,
			config=config
	}
    call nanopolish {
        input:
            reads=reads,
            allguppy=guppy.allguppy
    }
    call minimapalign {
        input:
            genome=genome,
            allguppy=guppy.allguppy
    }
    call filter {
        input:
            sortedbam=minimapalign.sortedbam
    }
    call methylation {
        input:
            filteredbam=filter.filteredbam,
            genome=genome,
            reads=reads,
            type=type,
            nanoindex=nanopolish.nanoindex,
            allguppy=guppy.allguppy,
            filteredbai=filter.filteredbai,
            nanofastaindex=nanopolish.nanofastaindex,
            nanoindexgzi=nanopolish.nanoindexgzi,
            nanoreaddb=nanopolish.nanoreaddb
    }
    call nanotobed {
        input:
            sortedbed=sortedbed,
            typeagain=typeagain,
            pythonscript=pythonscript,
            methylationcalls=methylation.methylationcalls
    }
    call bedtobigwig {
        input:
            chromsizes=chromsizes,
            FivemCavgbedgraph=nanotobed.FivemCavgbedgraph
    }
    output {
        File allguppy = guppy.allguppy
        File filteredbam = filter.filteredbam
        File methylationfreq = methylation.methylationcalls
        File FivemCbed = nanotobed.FivemCbed
        File FivemCavgbedgraph = nanotobed.FivemCavgbedgraph
        File FivemCavgbw = bedtobigwig.FivemCavgbw
	}

	meta {
		author: ""
		email:""
		description: "Workflow for using guppy and nanopolish to process nanopore output"
	}
}
task guppy {
    input {
        File reads
        File config
        File model
    }
    command <<<
    mkdir ./samples
    tar zxvf ~{reads} --strip-components=1 -C ./samples
    mkdir ./out
    guppy_basecaller -i ./samples -s ./out -c ~{config} -m ~{model} -x "cuda:all" --num_callers 4 --gpu_runners_per_device 8
    cat ./out/**/*.fastq > ./out/allguppy.fastq
    >>>
    runtime {
        gpuType: "nvidia-tesla-v100"
        gpuCount: 1
        nvidiaDriverVersion: "470.161.03"
        zones: ["us-central1-a"] 
		docker: "us-central1-docker.pkg.dev/aryeelab/docker/megalodon"
		memory: "64 GB"
		disks: "local-disk 1000 SSD"
		cpu: 12
    }
    output {
        File allguppy = "./out/allguppy.fastq"
    }
}

task nanopolish {
    input {
        File allguppy
        File reads
    }
    command <<<
    mkdir ./index
    cp ~{allguppy} ./index/
    mkdir ./samples
    tar zxvf ~{reads} --strip-components=1 -C ./samples
    nanopolish index -d ./samples ./index/allguppy.fastq
    >>>
    runtime {
        docker: "us-central1-docker.pkg.dev/aryeelab/docker/nanopolish:latest"
		memory: "64G"
		disks: "local-disk 1000 SSD"
		cpu: 16 
    }
    output {
        File nanoindex = "./index/allguppy.fastq.index"
        File nanofastaindex = "./index/allguppy.fastq.index.fai"
        File nanoindexgzi = "./index/allguppy.fastq.index.gzi"
        File nanoreaddb = "./index/allguppy.fastq.index.readdb"
    }
}
task minimapalign {
    input {
        File allguppy
        File genome
    }
    command <<<
    minimap2 -a -x map-ont ~{genome} ~{allguppy} | samtools sort -T tmp -o sorted.bam
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
task methylation {
    input {
        File reads
        File allguppy
        File filteredbam
        File filteredbai
        File nanoindex
        File nanofastaindex
        File nanoindexgzi
        File nanoreaddb
        File genome
        String type
    }
    command <<<
    mkdir ./samples
    tar zxvf ~{reads} --strip-components=1 -C ./samples
    mkdir ./temp
    cp ~{allguppy} ~{filteredbam} ~{filteredbai} ~{nanoindex} ~{nanofastaindex} ~{nanoindexgzi} ~{nanoreaddb} ./temp/
    nanopolish call-methylation --methylation=~{type} -t 8 -r ./temp/allguppy.fastq -b ./temp/filtered.bam -g ~{genome} > methylationcalls.tsv
    >>>
    runtime {
        docker: "us-central1-docker.pkg.dev/aryeelab/docker/nanopolish:latest"
		memory: "64G"
		disks: "local-disk 500 SSD"
		cpu: 8
    }
    output {
        File methylationcalls = "methylationcalls.tsv"
    }
}
task nanotobed {
    input {
        File pythonscript
        File methylationcalls
        File sortedbed
        String typeagain
    }
    command <<<
    python3 ~{pythonscript} -s -m ~{typeagain} ~{methylationcalls} > methylationfrequency.tsv
    tail -n +2 methylationfrequency.tsv | awk '{ print $1"\t"$2"\t"$3+1"\tid-"NR"\t"$7; }' | sort-bed - > FivemC.percentage.bed
    bedops --chop 1000 ~{sortedbed} | bedmap --faster --echo --mean --count --delim "\t" --skip-unmapped - FivemC.percentage.bed | cat | cut -f 1,2,3,4 | sort -k1,1 -k2,2n > nanopolish5mC.1k.bedgraph
    >>>
    runtime {
        docker: "us-central1-docker.pkg.dev/aryeelab/docker/bedops:latest"
		memory: "64G"
		disks: "local-disk 500 SSD"
		cpu: 8
    }
    output {
        File FivemCbed = "FivemC.percentage.bed"
        File FivemCavgbedgraph = "nanopolish5mC.1k.bedgraph"
    }
}
task bedtobigwig {
    input {
        File FivemCavgbedgraph
        File chromsizes
    }
    command <<<
    bedGraphToBigWig ~{FivemCavgbedgraph} ~{chromsizes} FivemCavg.bw
    >>>
    runtime {
        docker: "us-central1-docker.pkg.dev/aryeelab/docker/bedtools:latest"
		memory: "64G"
		disks: "local-disk 500 SSD"
		cpu: 8
    }
    output {
        File FivemCavgbw = "FivemCavg.bw"
    }
}