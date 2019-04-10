workflow preprocess_flowcell {

    String run_id
    String working_dir
    String image_dir
    File fast5_zip
    File ref_genome

    # Basecall and demultiplex with albacore
    call basecall_and_demultiplex {input: run_id = run_id, fast5_zip = fast5_zip, working_dir = working_dir,image_dir = image_dir}

    scatter (fastq_gz in basecall_and_demultiplex.fastq_gzs) {
    call removeReadsWithDuplicateID {input: fastq_gz = fastq_gz}
    # Align with minimap2
    call align {input: fastq_gz = removeReadsWithDuplicateID.dedup_fastq_gz, ref_genome = ref_genome, working_dir = working_dir,image_dir = image_dir}
    # Call methylation with Nanopolish
    call call_methylation {input:   fast5_zip = fast5_zip,
    					sequence_summary = basecall_and_demultiplex.sequence_summary,
                                        fastq_gz = removeReadsWithDuplicateID.dedup_fastq_gz,
                                        bam = align.bam,
                                        bai = align.bai,
                                        ref_genome = ref_genome,
					working_dir = working_dir,
					image_dir = image_dir}

     # Summarize methylation by read
     call methylation_by_read {input: base_methylation_calls = call_methylation.methylation_calls, working_dir = working_dir, image_dir = image_dir}

    }

    call demux_sample_sheet {input: fastq_gzs = basecall_and_demultiplex.fastq_gzs,
                                    bams = align.bam,
                                    bais = align.bai,
                                    methylation_calls = call_methylation.methylation_calls,
                                    read_methylation_calls = methylation_by_read.read_methylation_calls,
				    working_dir = working_dir,
				    image_dir = image_dir
                           }
    
}

task basecall_and_demultiplex {
	String run_id
	File fast5_zip
    	String flowcell_id
    	String kit_id
    	Int min_qscore
	String working_dir
	String image_dir
	String pattern = '${working_dir}/*.fq.gz'

	command <<<
		module load local/singularity/2.5.2
		fast5_path=`jar tf ${fast5_zip} | grep 'fast5/$'` # find path to fast5 dir within zip
	        echo "Unzipping ${fast5_zip}"
        	jar -xf ${fast5_zip}
        	echo "Moving $fast5_path to `basename $fast5_path`"
        	mv $fast5_path ${working_dir}/fast5
        	# Basecall
		singularity exec --bind ${working_dir}:${working_dir} ${image_dir}/guppy.simg bash -c "echo step1  && \
guppy_basecaller -r -i ${working_dir}/fast5 -s guppy_basecaller -q 0 --flowcell ${flowcell_id} --kit ${kit_id} --cpu_threads_per_caller 2 --qscore_filtering --min_qscore ${min_qscore} && \
                cat ${working_dir}/guppy_basecaller/guppy_basecaller_log* > guppy_basecaller.log && \
		guppy_barcoder -i ${working_dir}/guppy_basecaller/pass -s ${working_dir}/guppy_barcoder --barcode_kits ${kit_id}"


		barcodes="`cd ${working_dir}/guppy_barcoder && \ls -d barcode0*`"
                for barcode in $barcodes;
                do
                        fastq_gz=${run_id}__$barcode.fq.gz
                        echo $fastq_gz
                        cat ${working_dir}/guppy_barcoder/$barcode/*.fastq | gzip -c >  ${working_dir}/$fastq_gz
                done
	 >>>
    output {
        File sequence_summary = "${working_dir}/guppy_basecaller/sequencing_summary.txt"
        File guppy_basecaller_log = "${working_dir}/guppy_basecaller.log"
        Array[File] fastq_gzs = glob(pattern)
    }
	
}


task removeReadsWithDuplicateID {
    File fastq_gz
    String base = basename(fastq_gz, ".fq.gz")
    String fastq = "${base}.fq"
    String dedup_fastq = "${base}.dedup.fq"

    command <<<
		zcat ${fastq_gz} | \
		awk '{ if (NR % 4 == 0) {printf("%s\n",$0)} else {printf("%s\t",$0) } }' | \
		awk '!seen[$0]++' | \
		sed 's/\t/\r\n/g' | \
		gzip > ${base}.dedup.fq.gz
	>>>
    output {
        File dedup_fastq_gz = "${base}.dedup.fq.gz"
    }
}


task align {

    File fastq_gz
    File ref_genome
    String base = basename(fastq_gz, ".fq.gz")
    String working_dir
    String image_dir
    command <<<
	module load local/singularity/2.5.2
	singularity exec --bind ${working_dir}:${working_dir} ${image_dir}/nanopore_minimap2.simg bash -c "cd ${working_dir} && \
	minimap2 -ax map-ont -t 1 ${ref_genome} ${fastq_gz} | samtools sort -o ${working_dir}/${base}.bam && \
        samtools index ${working_dir}/${base}.bam "
    >>>
    
    output {
        File bam = "${working_dir}/${base}.bam"
        File bai = "${working_dir}/${base}.bam.bai"
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
    String working_dir
    String image_dir

    command <<<
	module load local/singularity/2.5.2

        # Unzip fast5
        fast5_path=`jar tf ${fast5_zip} | grep 'fast5/$'` # find path to fast5 dir within zip
        echo "Unzipping ${fast5_zip}"
        jar -xf ${fast5_zip}

	singularity exec --bind ${working_dir}:${working_dir} ${image_dir}/nanopolish-v0.11.0.simg bash -c "echo nanopolish && \
        /nanopolish/nanopolish index -d $fast5_path -s ${sequence_summary} ${fastq_gz} && \
        /nanopolish/nanopolish call-methylation -t 1 -r ${fastq_gz} -b ${bam} -g ${ref_genome} > ${working_dir}/${base}.methylation_calls.tsv"
#	/nanopolish/scripts/calculate_methylation_frequency.py -i ${working_dir}/${base}.methylation_calls.tsv > ${working_dir}/${base}.methylation_frequency.tsv"
    >>>

    output {
        File methylation_calls = "${working_dir}/${base}.methylation_calls.tsv"
#        File methylation_frequencies = "${working_dir}/${base}.methylation_frequency.tsv"
    }

}

task methylation_by_read {

    File base_methylation_calls
    String base = basename(base_methylation_calls, ".methylation_calls.tsv")
    String working_dir
    String image_dir

    command <<<
	module load local/singularity/2.5.2
        singularity exec --bind ${working_dir}:${working_dir} ${image_dir}/nanopore_util.simg bash -c " cp ${base_methylation_calls} . && \
        perl /usr/local/bin/methylation_by_read.pl `basename ${base_methylation_calls}` "
    >>>

    output {
        File read_methylation_calls = "${base}.read-methylation-calls.tsv"
    }

}

task demux_sample_sheet {
    Array[File] fastq_gzs
    Array[String] bams
    Array[String] bais
    Array[String] methylation_calls
    Array[String] read_methylation_calls
    String working_dir 
    String image_dir 

    command <<<
        echo fastq_gz ${sep=' ' fastq_gzs} >> samples_t.txt
        echo bam ${sep=' ' bams} >> samples_t.txt
        echo bai ${sep=' ' bais} >> samples_t.txt
        echo methylation_calls ${sep=' ' methylation_calls} >> samples_t.txt
        echo read_methylation_calls ${sep=' ' read_methylation_calls} >> samples_t.txt


        module load local/singularity/2.5.2
        singularity exec --bind ${working_dir}:${working_dir} ${image_dir}/nanopore_util.simg bash -c " cat samples_t.txt | rs -c' ' -C',' -T > samples.csv && \
        /usr/local/bin/add_flowcell_and_barcode_columns.R samples.csv samples.csv"
    >>>

    output {
        File samples = "samples.csv"
    }
}

