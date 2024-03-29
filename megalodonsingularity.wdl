version 1.0

workflow megalodon {

	input {
		String sampledir
		String device
		String modmotifs
		String outdir
		File sortedchromsize
		File genome
		File chromsizes
	}

	call meg {
		input:
			sampledir=sampledir,
			outdir=outdir,
			device=device,
			modmotifs=modmotifs,
			genome=genome
	}
	call QC {
		input:
			outbam=meg.outbam
	}
	# call bamtobed {
	# 	input:
	# 		outdir=outdir
	# }
	# call bedtobw {
	# 	input:
	# 		chromsizes=chromsizes
	# 		mappingsbed=bamtobed.mappingsbed
	# }
	call bedtobedgraph {
		input:
			outdir=outdir,
			FivemCbed=meg.FivemCbed,
			SixmAbed=meg.SixmAbed
	}
	call smoothing {
		input:
			sortedchromsize=sortedchromsize,
			FivemCpercentbed=bedtobedgraph.FivemCpercentbed,
			FivemCcoveragebed=bedtobedgraph.FivemCcoveragebed,
			SixmApercentbed=bedtobedgraph.SixmApercentbed,
			SixmAcoveragebed=bedtobedgraph.SixmAcoveragebed
	}
	call bedgraphtobigwig {
		input:
			outdir=outdir,
			FivemCpercentbedgraph=bedtobedgraph.FivemCpercentbedgraph,
			FivemCcoveragebedgraph=bedtobedgraph.FivemCcoveragebedgraph,
			SixmApercentbedgraph=bedtobedgraph.SixmApercentbedgraph,
			SixmAcoveragebedgraph=bedtobedgraph.SixmAcoveragebedgraph,
			chromsizes=chromsizes,
			FivemCpercentavgbedgraph=smoothing.FivemCpercentavgbedgraph,
			FivemCcoverageavgbedgraph=smoothing.FivemCcoverageavgbedgraph,
			SixmApercentavgbedgraph=smoothing.SixmApercentavgbedgraph,
			SixmAcoverageavgbedgraph=smoothing.SixmAcoverageavgbedgraph
	}
	output {
		File FivemCpercentbg = bedtobedgraph.FivemCpercentbedgraph
		File FivemCcoveragebg = bedtobedgraph.FivemCcoveragebedgraph
		File SixmApercentbg = bedtobedgraph.SixmApercentbedgraph
		File SixmAcoveragebg = bedtobedgraph.SixmAcoveragebedgraph
		File FivemCpercentbw = bedgraphtobigwig.FivemCpercentbw
		File FivemCcoveragebw = bedgraphtobigwig.FivemCcoveragebw
		File SixmApercentbw = bedgraphtobigwig.SixmApercentbw
		File SixmAcoveragebw = bedgraphtobigwig.SixmAcoveragebw
		File FivemCpercentsmooth = bedgraphtobigwig.FivemCpercentsmooth
		File FivemCcoveragesmooth = bedgraphtobigwig.FivemCcoveragesmooth
		File SixmApercentsmooth = bedgraphtobigwig.SixmApercentsmooth
		File SixmAcoveragesmooth = bedgraphtobigwig.SixmAcoveragesmooth
		File basecoverage = QC.base_coverage
		String reads = QC.num_reads
		# File mappingsbw = bedtobw.mappingsbw
	}

	meta {
		author: ""
		email:""
		description: "Workflow for turning nanopore sequencing output data into bigwigs"
	}
}
task meg {
	input {
		String sampledir
		String outdir
		File genome
		String modmotifs
		String device
	}
	command <<<
		mkdir ./out
		megalodon \
		~{sampledir} \
		--output-directory ~{outdir} \
		--overwrite \
		--guppy-server-path /usr/bin/guppy_basecall_server \
		--guppy-params "-d /aryeelab/nanopore/rerio/basecall_models/" \
		--guppy-config res_dna_r941_min_modbases-all-context_v001.cfg \
		--outputs basecalls mappings mod_mappings mods per_read_mods \
		--reference ~{genome} ~{modmotifs} \
		--devices ~{device} --mod-min-prob 0
	>>>
	runtime {
		image: "/aryeelab/singularity/megalodon.sif"
		bindpath: "/aryeelab"
		host:"/misc"
	}
	output {
		File FivemCbed = "${outdir}/modified_bases.5mC.bed"
		File SixmAbed = "${outdir}/modified_bases.6mA.bed"
		File outbam = "${outdir}/mappings.bam"
	}
}
task QC {
    input {
        File outbam
    }
    command <<<
        samtools view -c ~{outbam}
        samtools sort ~{outbam} -o bamoutputsorted.bam
        samtools index bamoutputsorted.bam
        echo "coverage number_of_bases" >> base_coverage.txt
        cov=$(samtools mpileup bamoutputsorted.bam | awk -v X="1" '$4>=X' | wc -l); echo 1 $cov >> base_coverage.txt
        cov=$(samtools mpileup bamoutputsorted.bam | awk -v X="2" '$4>=X' | wc -l); echo 2 $cov >> base_coverage.txt
        cov=$(samtools mpileup bamoutputsorted.bam | awk -v X="5" '$4>=X' | wc -l); echo 5 $cov >> base_coverage.txt
        cov=$(samtools mpileup bamoutputsorted.bam | awk -v X="10" '$4>=X' | wc -l); echo 10 $cov >> base_coverage.txt
        cov=$(samtools mpileup bamoutputsorted.bam | awk -v X="20" '$4>=X' | wc -l); echo 20 $cov >> base_coverage.txt
        cov=$(samtools mpileup bamoutputsorted.bam | awk -v X="30" '$4>=X' | wc -l); echo 30 $cov >> base_coverage.txt
    >>>
    runtime {
		image: "/misc/aryeelab-data/samtools.sif"
		bindpath: "/aryeelab"
		host:"/misc"
	}
    output {
        String num_reads = read_lines(stdout())[0]
        File base_coverage = "base_coverage.txt"
    }
}
# task bamtobed {
# 	input {
# 		String outdir
# 	}
# 	command <<<
# 	bedtools bamtobed -i ~{outdir}/mappings.bam | sort -k1,1 -k2,2n > mappings.bed
# 	>>>
# 	runtime {
# 		docker: ":"
# 	}
# 	output {
# 		File mappingsbed = "mappings.bed"
# 	}
# }
# task bedtobw {
# 	input {
# 		File mappingsbed
# 		File chromsizes
# 	}
# 	command <<<
# 	bedGraphToBigWig ~{mappingsbed} ~{chromsizes} mappings.bw
# 	>>>
# 	runtime {
# 		docker: ":"
# 	}
# 	output {
# 		File mappingsbw = "mappings.bw"
# 	}
# }
task bedtobedgraph {
	input {
		String outdir
		File FivemCbed
		File SixmAbed
	}
	command <<<
		cat ~{FivemCbed} | awk '{ print $1"\t"$2"\t"$3"\tid-"NR"\t"$11; }' | sort-bed - > FivemC.percentage.bed
		cat ~{FivemCbed} | awk '{ print $1"\t"$2"\t"$3"\tid-"NR"\t"$10; }' | sort-bed - > FivemC.coverage.bed
		cat ~{SixmAbed} | awk '{ print $1"\t"$2"\t"$3"\tid-"NR"\t"$11; }' | sort-bed - > SixmA.percentage.bed
		cat ~{SixmAbed} | awk '{ print $1"\t"$2"\t"$3"\tid-"NR"\t"$10; }' | sort-bed - > SixmA.coverage.bed
		cat ~{FivemCbed} | cut -f 1,2,3,11 | sort -k1,1 -k2,2n > FivemC.percent.bedgraph
		cat ~{FivemCbed} | cut -f 1,2,3,10 | sort -k1,1 -k2,2n > FivemC.coverage.bedgraph
		cat ~{SixmAbed} | cut -f 1,2,3,11 | sort -k1,1 -k2,2n > SixmA.percent.bedgraph
		cat ~{SixmAbed} | cut -f 1,2,3,10 | sort -k1,1 -k2,2n > SixmA.coverage.bedgraph
	>>>
	runtime {
		image: "/misc/aryeelab-data/bedops.sif"
		bindpath: "/aryeelab"
		host:"/misc"
	}
	output {
		File FivemCpercentbed = "FivemC.percentage.bed"
		File FivemCcoveragebed = "FivemC.coverage.bed"
		File SixmApercentbed = "SixmA.percentage.bed"
		File SixmAcoveragebed = "SixmA.coverage.bed"
		File FivemCpercentbedgraph = "FivemC.percent.bedgraph"
		File FivemCcoveragebedgraph = "FivemC.coverage.bedgraph"
		File SixmApercentbedgraph = "SixmA.percent.bedgraph"
		File SixmAcoveragebedgraph = "SixmA.coverage.bedgraph"
	}	
}
task smoothing {
	input {
		File sortedchromsize
		File FivemCpercentbed
		File FivemCcoveragebed
		File SixmApercentbed
		File SixmAcoveragebed
	}

	command <<<
		bedops --chop 1000 ~{sortedchromsize} | bedmap --faster --echo --mean --count --delim "\t" --skip-unmapped - ~{FivemCpercentbed} | cat | cut -f 1,2,3,4 | sort -k1,1 -k2,2n > FivemCpercent.avg.bedgraph
		bedops --chop 1000 ~{sortedchromsize} | bedmap --faster --echo --mean --count --delim "\t" --skip-unmapped - ~{FivemCcoveragebed} | cat | cut -f 1,2,3,4 | sort -k1,1 -k2,2n > FivemCcoverage.avg.bedgraph
		bedops --chop 1000 ~{sortedchromsize} | bedmap --faster --echo --mean --count --delim "\t" --skip-unmapped - ~{SixmApercentbed} | cat | cut -f 1,2,3,4 | sort -k1,1 -k2,2n > SixmApercent.avg.bedgraph
		bedops --chop 1000 ~{sortedchromsize} | bedmap --faster --echo --mean --count --delim "\t" --skip-unmapped - ~{SixmAcoveragebed} | cat | cut -f 1,2,3,4 | sort -k1,1 -k2,2n > SixmAcoverage.avg.bedgraph
	>>>

	runtime {
		image: "/misc/aryeelab-data/bedops.sif"
		bindpath: "/aryeelab"
		host:"/misc"
	}

	output {
		File FivemCpercentavgbedgraph = "FivemCpercent.avg.bedgraph"
		File FivemCcoverageavgbedgraph = "FivemCcoverage.avg.bedgraph"
		File SixmApercentavgbedgraph = "SixmApercent.avg.bedgraph"
		File SixmAcoverageavgbedgraph = "SixmAcoverage.avg.bedgraph"
	}
}
task bedgraphtobigwig {
	input {
		String outdir
		File FivemCpercentbedgraph
		File FivemCcoveragebedgraph
		File SixmApercentbedgraph
		File SixmAcoveragebedgraph
		File FivemCpercentavgbedgraph
		File FivemCcoverageavgbedgraph
		File SixmApercentavgbedgraph
		File SixmAcoverageavgbedgraph
		File chromsizes
	}
	command <<<
		bedGraphToBigWig ~{FivemCpercentbedgraph} ~{chromsizes} FivemC.percentage.bw
		bedGraphToBigWig ~{SixmApercentbedgraph} ~{chromsizes} SixmA.percentage.bw
		bedGraphToBigWig ~{FivemCcoveragebedgraph} ~{chromsizes} FivemC.coverage.bw
		bedGraphToBigWig ~{SixmAcoveragebedgraph} ~{chromsizes} SixmA.coverage.bw
		bedGraphToBigWig ~{FivemCpercentavgbedgraph} ~{chromsizes} FivemC.percentage.smoothed.bw
		bedGraphToBigWig ~{FivemCcoverageavgbedgraph} ~{chromsizes} FivemC.coverage.smoothed.bw
		bedGraphToBigWig ~{SixmApercentavgbedgraph} ~{chromsizes} SixmA.percentage.smoothed.bw
		bedGraphToBigWig ~{SixmAcoverageavgbedgraph} ~{chromsizes} SixmA.coverage.smoothed.bw
	>>>
	runtime {
		image: "/misc/aryeelab-data/bedgraphtobigwig.sif"
		bindpath: "/aryeelab"
		host:"/misc"
	}
	output {
		File FivemCpercentbw = "FivemC.percentage.bw"
		File FivemCcoveragebw = "FivemC.coverage.bw"
		File SixmApercentbw = "SixmA.percentage.bw"
		File SixmAcoveragebw = "SixmA.coverage.bw"
		File FivemCpercentsmooth = "FivemC.percentage.smoothed.bw"
		File FivemCcoveragesmooth = "FivemC.coverage.smoothed.bw"
		File SixmApercentsmooth = "SixmA.percentage.smoothed.bw"
		File SixmAcoveragesmooth = "SixmA.coverage.smoothed.bw"
	}
}