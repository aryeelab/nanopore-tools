workflow preprocess_flowcell {

    String version = "dev"
    #File monitoring_script = "gs://aryeelab/scripts/monitor_v2.sh"
    File monitoring_script = "monitor_v2.sh"
    String flowcell_id
    File fast5_zip
    File ref_genome
    File cpg_islands
    File chrs
    File compartments

    Int min_reads_per_barcode = 100

    Int disk_size = ceil(size(fast5_zip, "GB")) * 10 + 20


    # Basecall and demultiplex with albacore
    call basecall_and_demultiplex {input:   flowcell_id = flowcell_id, 
                                            fast5_zip = fast5_zip, 
                                            min_reads_per_barcode = min_reads_per_barcode, 
                                            disk_size = disk_size,
                                            monitoring_script = monitoring_script}

    scatter (fastq_gz in basecall_and_demultiplex.fastq_gzs) {
        call removeReadsWithDuplicateID {input: fastq_gz = fastq_gz}
        
        # Align with minimap2
        call align {input: fastq_gz = removeReadsWithDuplicateID.dedup_fastq_gz, 
                           ref_genome = ref_genome, 
                           monitoring_script = monitoring_script}
        
        # Call methylation with Nanopolish
        call call_methylation {input:   fast5_zip = fast5_zip,
                                        sequence_summary = basecall_and_demultiplex.sequence_summary,
                                        fastq_gz = removeReadsWithDuplicateID.dedup_fastq_gz,
                                        bam = align.bam,
                                        bai = align.bai,
                                        ref_genome = ref_genome,
                                        disk_size = disk_size,
                                        monitoring_script = monitoring_script}

         # Summarize methylation by read
         call methylation_by_read {input: base_methylation_calls = call_methylation.methylation_calls, version = version, threshold_ll = 2.5, cpg_islands = cpg_islands, chrs = chrs, compartments = compartments}
    }

    call demux_sample_sheet {input: flowcell_id = flowcell_id,
                                    fastq_gzs = basecall_and_demultiplex.fastq_gzs,
                                    dedup_fastq_gzs = removeReadsWithDuplicateID.dedup_fastq_gz,
                                    bams = align.bam,
                                    bais = align.bai,
                                    methylation_calls = call_methylation.methylation_calls,
                                    read_methylation_calls = get_methylation_by_read.read_methylation_calls,
                                   # methylation_read_sequence = get_methylation_read_sequence.reads,
                                    version = version}
    
    output {
        File guppy_basecaller_log = basecall_and_demultiplex.guppy_basecaller_log
        File sequence_summary = basecall_and_demultiplex.sequence_summary
        File sample_sheet = demux_sample_sheet.samples
        String pipeline_version = "${version}"
    }

}

task basecall_and_demultiplex {
	String flowcell_id
	File fast5_zip
    String flowcell_type_id
    String kit_id
    String device = "auto"
    Int min_qscore
    Int min_reads_per_barcode
    Boolean gpu = true
    Int disk_size
    File monitoring_script

	command <<<
    	chmod u+x ${monitoring_script}
        ${monitoring_script} > monitoring.log &
	
	    extension=`echo ${fast5_zip} | sed 's/.*\.//'`
	    if [ $extension == "zip" ]
	    then    
            fast5_path=`jar tf ${fast5_zip} | grep 'fast5/$'` # find path to fast5 dir within zip
            echo "Unzipping ${fast5_zip}"
            jar -xf ${fast5_zip}
            echo "Moving $fast5_path to `basename $fast5_path`"
            mv $fast5_path fast5
	    elif [ $extension == "tar" ]
	    then    
            echo "Untarring ${fast5_zip}"
            mkdir fast5
            tar -xf ${fast5_zip} -C fast5
        else 
            echo "ERROR: Input file does not end in tar or zip."
            exit 1
        fi
        
        # Basecall
		guppy_basecaller -r -i fast5 -s guppy_basecaller -q 0 --flowcell ${flowcell_type_id} --kit ${kit_id} --qscore_filtering --min_qscore ${min_qscore} ${if gpu then '--device' else ''} ${if gpu then device else ''} 
		    
        cat guppy_basecaller/guppy_basecaller_log* > guppy_basecaller.log
		guppy_barcoder -i guppy_basecaller/pass -s guppy_barcoder
		\ls -d guppy_barcoder/*
		barcodes="`cd guppy_barcoder && \ls -d barcode0* unclassified`"
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
        continueOnReturnCode: false
        docker: if gpu then "quay.io/aryeelab/guppy-gpu" else "quay.io/aryeelab/guppy-cpu"
        bootDiskSizeGb: 200
        disks: "local-disk ${disk_size} HDD"
        gpuType: "nvidia-tesla-p100"
        gpuCount: 1
        zones: "us-central1-c"
        preemptible: 1
        simg: "${if gpu then 'guppy-gpu.simg' else 'guppy-cpu.simg'}"
    }
    output {
        File sequence_summary = "guppy_basecaller/sequencing_summary.txt"
        File guppy_basecaller_log = "guppy_basecaller.log"
        Array[File] fastq_gzs = glob("*.fq.gz")
        File monitoring_log = "monitoring.log"
    }
	
}


task removeReadsWithDuplicateID {
    File fastq_gz
    String base = basename(fastq_gz, ".fq.gz")
    String fastq = "${base}.fq"
    String dedup_fastq = "${base}.dedup.fq"
    Int disk_size = ceil(size(fastq_gz, "GB")) * 3 + 20

    command <<<
		zcat ${fastq_gz} | \
		awk '{ if (NR % 4 == 0) {printf("%s\n",$0)} else {printf("%s\t",$0) } }' | \
		awk '!seen[$0]++' | \
		sed 's/\t/\n/g' | \
		gzip > ${base}.dedup.fq.gz
	>>>
	
    runtime {
        continueOnReturnCode: false
        docker: "quay.io/aryeelab/nanopore-util"
        memory: "30GB"
        disks: "local-disk ${disk_size} HDD"
        simg: "nanopore-util.simg"
    }
    
    output {
        File dedup_fastq_gz = "${base}.dedup.fq.gz"
    }
}


task align {
    File fastq_gz
    File ref_genome
    String base = basename(fastq_gz, ".fq.gz")
    Int disk_size = ceil(size(fastq_gz, "GB")) * 3 + 20
    File monitoring_script

    command <<<
        chmod u+x ${monitoring_script}
        ${monitoring_script} > monitoring.log &
        
        minimap2 -ax map-ont -t 1 ${ref_genome} ${fastq_gz} | samtools view -b -F 2304 | samtools sort -o ${base}.bam
        samtools index ${base}.bam    
    >>>
    
    runtime {
        continueOnReturnCode: false
        docker: "quay.io/aryeelab/minimap2"
        memory: "100GB"
        disks: "local-disk ${disk_size} HDD"
        simg: "minimap2.simg"
    } 
       
    output {
        File bam = "${base}.bam"
        File bai = "${base}.bam.bai"
        File monitoring_log = "monitoring.log"
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
    Int disk_size 
    File monitoring_script

    command <<<
        chmod u+x ${monitoring_script}
        ${monitoring_script} > monitoring.log &
    
        # Unzip fast5
        fast5_path=`jar tf ${fast5_zip} | grep 'fast5/$'` # find path to fast5 dir within zip
        echo "Unzipping ${fast5_zip}"
        jar -xf ${fast5_zip}

        /nanopolish/nanopolish index -d $fast5_path -s ${sequence_summary} ${fastq_gz}
        /nanopolish/nanopolish call-methylation -t 8 -r ${fastq_gz} -b ${bam} -g ${ref_genome} > ${base}.methylation_calls.tsv
    >>>
    
    runtime {
        continueOnReturnCode: false
        docker: "quay.io/aryeelab/nanopolish"
        cpu: 8
        disks: "local-disk ${disk_size} HDD"
        preemptible: 1
        simg: "nanopolish.simg"
    }
    
    output {
        File methylation_calls = "${base}.methylation_calls.tsv"
        File monitoring_log = "monitoring.log"
    }

}

task get_methylation_by_read {
    File base_methylation_calls
    File cpg_islands
    File chrs
    File compartments
    Int threshold_ll
    String base = basename(base_methylation_calls, ".methylation_calls.tsv")
    Int disk_size = ceil(size(base_methylation_calls, "GB")) * 2 + 20
    File monitoring_script

    command <<<
        chmod u+x ${monitoring_script}
        ${monitoring_script} > monitoring.log &
    	Rscript /usr/local/bin/kmer_to_read.R ${base_methylation_calls} ${base}.read-methylation-calls.tsv ${base}.annotated_kmers.tsv ${threshold_ll} ${cpg_islands} ${chrs} ${compartments}
    >>>
    
   runtime {
        continueOnReturnCode: false
        docker: "quay.io/aryeelab/nanopore-util"
        memory: "30GB"
        disks: "local-disk ${disk_size} HDD"
        simg: "nanopore-util.simg"
    }
    
   output {
        File read_methylation_calls = "${base}.read-methylation-calls.tsv"
        File monitoring_log = "monitoring.log"
        File annotated_kmers = "${base}.annotated_kmers.tsv"
    }

}

task get_methylation_read_sequence {
    File base_methylation_calls
    String base = basename(base_methylation_calls, ".methylation_calls.tsv")
    Int disk_size = ceil(size(base_methylation_calls, "GB")) * 2 + 20
    File monitoring_script

    command <<<
        chmod u+x ${monitoring_script}
        ${monitoring_script} > monitoring.log &

        python3.7 /usr/local/bin/methylation_read_sequence.py ${base_methylation_calls} ${base}.methylation_reads.tsv
    >>>

   runtime {
        continueOnReturnCode: false
        docker: "quay.io/aryeelab/nanopore-read-tools"
        disks: "local-disk ${disk_size} HDD"
        simg: "nanopore-read-tools.simg"
    }
    
   output {
        File reads = "${base}.methylation_reads.tsv"
        File monitoring_log = "monitoring.log"
    }
}


task demux_sample_sheet {
    String version
    String flowcell_id 
    Array[String] fastq_gzs
    Array[String] dedup_fastq_gzs
    Array[String] bams
    Array[String] bais
    Array[String] methylation_calls
    Array[String] read_methylation_calls

    command <<<
        echo fastq_gz ${sep=' ' fastq_gzs} >> samples_t.txt
        echo dedup_fastq_gz ${sep=' ' dedup_fastq_gzs} >> samples_t.txt
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
        continueOnReturnCode: false
        docker: "quay.io/aryeelab/nanopore-util"
        simg: "nanopore-util.simg"
    }
    
    output {
        File samples = "${flowcell_id}.samples.csv"
    }
}
