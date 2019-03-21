workflow preprocess_flowcell {

    String run_id
    File fast5_zip
    File ref_genome

    # Basecall and demultiplex with albacore
    call basecall_and_demultiplex {input: run_id = run_id, fast5_zip = fast5_zip}
    
    scatter (fastq_gz in basecall_and_demultiplex.fastq_gzs) {    
        # Remove reads with duplicate IDs
        call deduplicate {input: fastq_gz = fastq_gz}

        # Align with minimap2
        call align {input: fastq_gz = deduplicate.dedup_fastq_gz, ref_genome = ref_genome}
        
        # Call methylation with Nanopolish
        call call_methylation {input:   fast5_zip = fast5_zip, 
                                        sequence_summary = basecall_and_demultiplex.sequence_summary,
                                        fastq_gz = deduplicate.dedup_fastq_gz,
                                        bam = align.bam,
                                        bai = align.bai,
                                        ref_genome = ref_genome}
                              }

}

task basecall_and_demultiplex {
    String run_id
    File fast5_zip
    String flowcell_id
    String kit_id
    Int min_qscore
    
    command <<<
        # Unzip fast5
        fast5_path=`jar tf ${fast5_zip} | grep 'fast5/$'` # find path to fast5 dir within zip
        echo "Unzipping ${fast5_zip}"
        jar -xf ${fast5_zip}
        echo "Moving $fast5_path to `basename $fast5_path`"
        mv $fast5_path /

        # Basecall
        guppy_basecaller -r -i /fast5 -s guppy_basecaller -q 0 --flowcell ${flowcell_id} --kit ${kit_id} --cpu_threads_per_caller 2 --qscore_filtering --min_qscore ${min_qscore}
        cat guppy_basecaller/guppy_basecaller_log* > guppy_basecaller.log
        
        # Demultiplex
        guppy_barcoder -i guppy_basecaller/pass -s guppy_barcoder --barcode_kits ${kit_id}

        barcodes="`cd guppy_barcoder && ls -d barcode0*`"
        for barcode in $barcodes;
        do
            fastq_gz=${run_id}__$barcode.fq.gz
            echo $fastq_gz
            cat guppy_barcoder/$barcode/*.fastq | gzip -c > $fastq_gz
        done
        
    >>>

    runtime {
        docker: "aryeelab/guppy"
    }
    
    output {
        File sequence_summary = "guppy_basecaller/sequencing_summary.txt"
        File guppy_basecaller_log = "guppy_basecaller.log"
        Array[File] fastq_gzs = glob("*.fq.gz")
    }
}

task deduplicate {
    File fastq_gz
    String base = basename(fastq_gz, ".fq.gz")
    String fastq = "${base}.fq"
    String dedup_fastq = "${base}.dedup.fq"
     
    command <<<        
        gunzip -c ${fastq_gz} > ${fastq}
        
        # Dummy dedup step.
        cp ${fastq} ${dedup_fastq}
        
        gzip ${dedup_fastq}
    >>>

    runtime {
        docker: "debian:stretch"
    }
    
    output {
        File dedup_fastq_gz = "${base}.dedup.fq.gz"
    }       
}

task align {

    File fastq_gz
    File ref_genome
    String base = basename(fastq_gz, ".fq.gz")

    command <<<
        minimap2 -ax map-ont -t 1 "${ref_genome}" "${fastq_gz}" | samtools sort -o "${base}.bam";
        samtools index "${base}.bam"
    >>> 

    runtime {
        docker: "aryeelab/nanopore_minimap2"
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

    command <<<
        echo "hello"
        
        # Unzip fast5
        fast5_path=`jar tf ${fast5_zip} | grep 'fast5/$'` # find path to fast5 dir within zip
        echo "Unzipping ${fast5_zip}"
        jar -xf ${fast5_zip}
        echo "Moving $fast5_path to /fast5"
        mv $fast5_path /
        
        # Index the reads
        /nanopolish/nanopolish index -d /fast5 -s ${sequence_summary} ${fastq_gz}

        # Call methylation
        /nanopolish/nanopolish call-methylation -t 1 -r ${fastq_gz} -b ${bam} -g ${ref_genome} > ${base}.methylation_calls.tsv 
    >>>
    
    runtime {
        docker: "aryeelab/nanopolish:v0.11.0"
    }
    
    output {
        File methylation_calls = "${base}.methylation_calls.tsv"
    } 
    
}