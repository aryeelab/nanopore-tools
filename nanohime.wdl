version 1.0

workflow nanohime {

	input {
        File reads
        File genome
        File config
        File model
        File nanohimetar
        File pythonscript
        File sortedbed
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
    call eventalign {
        input:
        allguppy=guppy.allguppy,
        filteredbam=filter.filteredbam,
        filteredbai=filter.filteredbai,
        nanoindex=nanopolish.nanoindex,
        nanofastaindex=nanopolish.nanofastaindex,
        nanoindexgzi=nanopolish.nanoindexgzi,
        nanoreaddb=nanopolish.nanoreaddb,
        genome=genome
    }
    call hime {
        input:
            genome=genome,
            eventaligntxt=eventalign.eventaligntxt,
            nanohimetar=nanohimetar
    }
    call tobed {
        input:
            FivemCmethylation=hime.FivemCmethylation,
            SixmAmethylation=hime.SixmAmethylation,
            sortedbed=sortedbed,
            pythonscript=pythonscript
    }
    call bedtobigwig {
        input:
            FivemCavgbedgraph=tobed.FivemCavgbedgraph,
            SixmAavgbedgraph=tobed.SixmAavgbedgraph,
            chromsizes=chromsizes
    }
    output {
        File allguppy = guppy.allguppy
        File filteredbam = filter.filteredbam
        File FivemCbed = tobed.FivemCbed
        File SixmAbed = tobed.SixmAbed
        File FivemCavgbw = bedtobigwig.FivemCavgbw
        File SixmAavgbw = bedtobigwig.SixmAavgbw
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
    tar zxvf ~{reads} -C ./samples
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
    tar zxvf ~{reads} -C ./samples
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
task eventalign {
    input {
        File allguppy
        File filteredbam
        File filteredbai
        File nanoindex
        File nanofastaindex
        File nanoindexgzi
        File nanoreaddb
        File genome
    }
    command <<<
    mkdir ./temp
    cp ~{allguppy} ~{filteredbam} ~{filteredbai} ~{nanoindex} ~{nanofastaindex} ~{nanoindexgzi} ~{nanoreaddb} ./temp/
    nanopolish eventalign -n -t 16 --reads ./temp/allguppy.fastq --bam ./temp/filtered.bam --genome ~{genome} --scale-events > out.eventalign.txt
    >>>
    runtime {
        docker: "us-central1-docker.pkg.dev/aryeelab/docker/nanopolish:latest"
		memory: "64G"
		disks: "local-disk 500 SSD"
		cpu: 8
    }
    output {
        File eventaligntxt = "out.eventalign.txt"
    }
}
task hime {
    input {
        File genome
        File eventaligntxt
        File nanohimetar
    }
    command <<<
    mkdir ./temp
    tar zxvf ~{nanohimetar} -C ./temp
    perl ./temp/nanoHiMe-main/perl_script/upper.pl ~{genome} > ./temp/ref_upper.fa
    samtools faidx ./temp/ref_upper.fa
    ./temp/nanoHiMe-main/nanoHiMe 6mA ref_upper.fa ~{eventaligntxt} output.6mA 50 50
    echo -e "chromosome\tstrand\tstart\tend\tread\tlog_lik_ratio\tsequence" | cat - output.6mA.methylation.txt > SixmAmethylation.tsv
    ./temp/nanoHiMe-main/nanoHiMe mCG ref_upper.fa ~{eventaligntxt} output.mCG
    echo -e "chromosome\tstrand\tstart\tend\tread\tlog_lik_ratio\tsequence" | cat - output.mcG.methylation.txt > FivemCmethylation.tsv
    >>>
    runtime {
        docker: "us-central1-docker.pkg.dev/aryeelab/docker/nanohime:latest"
		memory: "64G"
		disks: "local-disk 500 SSD"
		cpu: 8
    }
    output {
        File SixmAmethylation = "SixmAmethylation.tsv"
        File FivemCmethylation = "FivemCmethylation.tsv"
    }
}
task tobed {
    input {
        File SixmAmethylation
        File FivemCmethylation
        File pythonscript
        File sortedbed
    }
    command <<<
    python3 ~{pythonscript} -s -m CG ~{FivemCmethylation} > FivemCmethylationfrequency.tsv
    tail -n +2 FivemCmethylationfrequency.tsv | awk '{ print $1"\t"$2"\t"$3+1"\tid-"NR"\t"$7; }' | sort-bed - > FivemC.percentage.bed
    bedops --chop 1000 ~{sortedbed} | bedmap --faster --echo --mean --count --delim "\t" --skip-unmapped - FivemC.percentage.bed | cat | cut -f 1,2,3,4 | sort -k1,1 -k2,2n > 5mC.1k.bedgraph
    python3 ~{pythonscript} -s -m A ~{SixmAmethylation} > SixmAmethylationfrequency.tsv
    tail -n +2 SixmAmethylationfrequency.tsv | awk '{ print $1"\t"$2"\t"$3+1"\tid-"NR"\t"$7; }' | sort-bed - > SixmA.percentage.bed
    bedops --chop 1000 ~{sortedbed} | bedmap --faster --echo --mean --count --delim "\t" --skip-unmapped - SixmA.percentage.bed | cat | cut -f 1,2,3,4 | sort -k1,1 -k2,2n > 6mA.1k.bedgraph
    >>>
    runtime {
        docker: "us-central1-docker.pkg.dev/aryeelab/docker/bedops:latest"
		memory: "64G"
		disks: "local-disk 500 SSD"
		cpu: 8
    }
    output {
        File FivemCbed = "FivemC.percentage.bed"
        File FivemCavgbedgraph = "5mC.1k.bedgraph"
        File SixmAbed = "SixmA.percentage.bed"
        File SixmAavgbedgraph = "6mA.1k.bedgraph"
    }
}
task bedtobigwig {
    input {
        File FivemCavgbedgraph
        File SixmAavgbedgraph
        File chromsizes
    }
    command <<<
    bedGraphToBigWig ~{FivemCavgbedgraph} ~{chromsizes} FivemCavg.bw
    bedGraphToBigWig ~{SixmAavgbedgraph} ~{chromsizes} SixmAavg.bw
    >>>
    runtime {
        docker: "us-central1-docker.pkg.dev/aryeelab/docker/bedtools:latest"
		memory: "64G"
		disks: "local-disk 500 SSD"
		cpu: 8
    }
    output {
        File FivemCavgbw = "FivemCavg.bw"
        File SixmAavgbw = "SixmAavg.bw"
    }
}