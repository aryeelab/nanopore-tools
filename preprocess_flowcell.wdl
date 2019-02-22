workflow preprocess_flowcell {

    # Run albacore to basecall and demultiplex
    call albacore

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
        Array[File] fastqs = glob("*.fq.gz")
    }
}