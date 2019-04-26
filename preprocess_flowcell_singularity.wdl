workflow preprocess_flowcell {

    String flowcell_id
    String image_dir
    File fast5_zip
    File ref_genome
    Int min_reads_per_barcode = 100

    # Basecall and demultiplex with albacore
    call basecall_and_demultiplex {input: flowcell_id = flowcell_id, fast5_zip = fast5_zip, min_reads_per_barcode = min_reads_per_barcode, image_dir = image_dir}

    scatter (fastq_gz in basecall_and_demultiplex.fastq_gzs) {
    call removeReadsWithDuplicateID {input: fastq_gz = fastq_gz, image_dir = image_dir}
    # Align with minimap2
    call align {input: fastq_gz = removeReadsWithDuplicateID.dedup_fastq_gz, ref_genome = ref_genome, image_dir = image_dir}
    # Call methylation with Nanopolish
    call call_methylation {input:   fast5_zip = fast5_zip,
    					sequence_summary = basecall_and_demultiplex.sequence_summary,
                                        fastq_gz = removeReadsWithDuplicateID.dedup_fastq_gz,
                                        bam = align.bam,
                                        bai = align.bai,
                                        ref_genome = ref_genome,
					image_dir = image_dir}

     # Summarize methylation by read
     call methylation_by_read {input: base_methylation_calls = call_methylation.methylation_calls, image_dir = image_dir}

    }

    call demux_sample_sheet {input: flowcell_id = flowcell_id,
                                    fastq_gzs = basecall_and_demultiplex.fastq_gzs,
                                    bams = align.bam,
                                    bais = align.bai,
                                    methylation_calls = call_methylation.methylation_calls,
                                    read_methylation_calls = methylation_by_read.read_methylation_calls,
				    image_dir = image_dir
                           }
    
}

task basecall_and_demultiplex {
	String flowcell_id
	File fast5_zip
    String flowcell_id
    String kit_id
    String device = "auto"
    Int min_qscore
    Int min_reads_per_barcode
	String image_dir

	command <<<
		fast5_path=`jar tf ${fast5_zip} | grep 'fast5/$'` # find path to fast5 dir within zip
	        echo "Unzipping ${fast5_zip}"
        	jar -xf ${fast5_zip}
        	echo "Moving $fast5_path to `basename $fast5_path`"
        	mv $fast5_path fast5

        	# Basecall
		guppy_basecaller -r -i fast5 -s guppy_basecaller -q 0 --flowcell ${flowcell_id} --kit ${kit_id} --device ${device} --qscore_filtering --min_qscore ${min_qscore}
                cat guppy_basecaller/guppy_basecaller_log* > guppy_basecaller.log
		guppy_barcoder -i guppy_basecaller/pass -s guppy_barcoder


		barcodes="`cd guppy_barcoder && \ls -d barcode0*`"
        	for barcode in $barcodes;
        	do
                    fastq_gz=${flowcell_id}__$barcode.fq.gz
                    echo "Processing $barcode"
                    numlines=$(cat guppy_barcoder/$barcode/*.fastq | wc -l)
                    numreads=$((numlines/4))
                    if [[ $numreads  -ge ${min_reads_per_barcode} ]]
                    then
                        cat guppy_barcoder/$barcode/*.fastq | gzip -c >  $fastq_gz
                        echo "Wrote $numreads reads to $fastq_gz"
                    fi
                    if [[ $numreads -lt ${min_reads_per_barcode} ]]
                        then
                        echo "Skipping since there are only $numreads reads (i.e. less than the ${min_reads_per_barcode}) read cutoff"
                    fi
                done

                # Exit with error if no barcodes are detected
                num_barcodes=$(echo -n $barcodes | wc -m)
                if [ $num_barcodes -eq 0 ]
                then
                    exit 1
                fi               
	 >>>
    runtime {
        backend: "singularity"
        continueOnReturnCode: false
        simg: "${image_dir}/guppy-gpu.simg"
    }
    output {
        File sequence_summary = "guppy_basecaller/sequencing_summary.txt"
        File guppy_basecaller_log = "guppy_basecaller.log"
        Array[File] fastq_gzs = glob("*.fq.gz")
    }
	
}


task removeReadsWithDuplicateID {
    File fastq_gz
    String base = basename(fastq_gz, ".fq.gz")
    String fastq = "${base}.fq"
    String dedup_fastq = "${base}.dedup.fq"
    String image_dir

    command <<<
		zcat ${fastq_gz} | \
		awk '{ if (NR % 4 == 0) {printf("%s\n",$0)} else {printf("%s\t",$0) } }' | \
		awk '!seen[$0]++' | \
		sed 's/\t/\r\n/g' | \
		gzip > ${base}.dedup.fq.gz
	>>>
    runtime {
        backend: "singularity"
        continueOnReturnCode: false
        simg: "${image_dir}/nanopore_util.simg"
    }
    output {
        File dedup_fastq_gz = "${base}.dedup.fq.gz"
    }
}


task align {

    File fastq_gz
    File ref_genome
    String base = basename(fastq_gz, ".fq.gz")
    String image_dir
    command <<<
	minimap2 -ax map-ont -t 1 ${ref_genome} ${fastq_gz} | samtools sort -o ${base}.bam
        samtools index ${base}.bam
    >>>
    runtime {
        backend: "singularity"
        continueOnReturnCode: false
        simg: "${image_dir}/nanopore_minimap2.simg"
    }    
    output {
        File bam = "${base}.bam"
        File bai = "${base}.bam.bai"
    }
}

task call_methylation {

    File fast5_zip
    File sequence_summary
    File fastq_gz
    File bam
    File bai
    File ref_genome
    String base = basename(fastq_gz, ".fq.gz")
    String image_dir

    command <<<
        # Unzip fast5
        fast5_path=`jar tf ${fast5_zip} | grep 'fast5/$'` # find path to fast5 dir within zip
        echo "Unzipping ${fast5_zip}"
        jar -xf ${fast5_zip}

        /nanopolish/nanopolish index -d $fast5_path -s ${sequence_summary} ${fastq_gz}
        /nanopolish/nanopolish call-methylation -t 1 -r ${fastq_gz} -b ${bam} -g ${ref_genome} > ${base}.methylation_calls.tsv
    >>>
    runtime {
        backend: "singularity"
        continueOnReturnCode: false
        simg: "${image_dir}/nanopolish.simg"
    }
    output {
        File methylation_calls = "${base}.methylation_calls.tsv"
    }

}

task methylation_by_read {

    File base_methylation_calls
    String base = basename(base_methylation_calls, ".methylation_calls.tsv")
    String image_dir

    command <<<
	Rscript /usr/local/bin/methylation_by_read.R ${base_methylation_calls} ${base}.read-methylation-calls.tsv
    >>>
   runtime {
        backend: "singularity"
        continueOnReturnCode: false
        simg: "${image_dir}/nanopore_util.simg"
    }
   output {
        File read_methylation_calls = "${base}.read-methylation-calls.tsv"
    }

}

task demux_sample_sheet {
    String flowcell_id 
    Array[File] fastq_gzs
    Array[String] bams
    Array[String] bais
    Array[String] methylation_calls
    Array[String] read_methylation_calls
    String image_dir 

    command <<<
        echo fastq_gz ${sep=' ' fastq_gzs} >> samples_t.txt
        echo bam ${sep=' ' bams} >> samples_t.txt
        echo bai ${sep=' ' bais} >> samples_t.txt
        echo methylation_calls ${sep=' ' methylation_calls} >> samples_t.txt
        echo read_methylation_calls ${sep=' ' read_methylation_calls} >> samples_t.txt


        cat samples_t.txt | datamash --output-delimiter=',' -t ' ' transpose > samples.csv 
        /usr/local/bin/add_flowcell_and_barcode_columns.R samples.csv samples.csv
        cp samples.csv "${flowcell_id}.samples.csv"
    	#cp samples.csv /data/aryee/sowmya/ctc_nanopore
	    #cp ${sep=' ' read_methylation_calls} /data/aryee/sowmya/ctc_nanopore/final_output_dir/
    >>>
    runtime {
        backend: "singularity"
        continueOnReturnCode: false
        simg: "${image_dir}/nanopore_util.simg"
    }
    output {
        File samples = "${flowcell_id}.samples.csv"
    }
}

