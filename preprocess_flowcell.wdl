workflow preprocess_flowcell {

    # Run albacore to basecall and demultiplex
    call albacore
    
    scatter (fastq_gz in albacore.fastq_gzs) {    
        # Remove reads with duplicate IDs
        call deduplicate {input: fastq_gz = fastq_gz}

        # Align with minimap2
        call align {input: fastq_gz = deduplicate.dedup_fastq_gz}
    }

}

task albacore {
    String run_id
    File fast5_zip
    String flowcell_id
    String kit_id
    
    command <<<
        # Unzip fast5
        fast5_path=`jar tf ${fast5_zip} | grep 'fast5/$'` # find path to fast5 dir within zip
        echo "Unzipping ${fast5_zip}"
        jar -xf ${fast5_zip}
        echo "Moving $fast5_path to `basename $fast5_path`"
        mv $fast5_path .

        read_fast5_basecaller.py -f  ${flowcell_id} -k ${kit_id} -t 1 -r -i  fast5 -o fastq -q 0 -s albacore

        barcodes=`ls albacore/workspace/pass`
        for barcode in $barcodes;
        do
            fastq_gz=${run_id}__$barcode.fq.gz
            echo $fastq_gz
            cat albacore/workspace/pass/$barcode/*.fastq | gzip -c > $fastq_gz
        done
        
    >>>

    runtime {
        docker: "aryeelab/nanopore_albacore"
    }
    
    output {
        File sequence_summary = "albacore/sequencing_summary.txt"
        File configuration = "albacore/configuration.cfg"
        File pipeline_log = "albacore/pipeline.log"
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
    File genome_index
    String base = basename(fastq_gz, ".fq.gz")

    command <<<
        minimap2 -ax map-ont -t 1 "${genome_index}" "${fastq_gz}" | samtools sort -o "${base}.bam";
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