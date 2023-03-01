version 1.0

workflow dorado_duplex {

    input {
        String sample_id
        File fast5_archive
        String basecall_model
    }

    call basecall_duplex {
        input:
            sample_id = sample_id,
            fast5_archive = fast5_archive,
            basecall_model = basecall_model
    }

    output {
        File duplex_bam = basecall_duplex.duplex_bam
        File pairs = basecall_duplex.pairs
    }

    meta {
        author: "Martin Aryee"
        email:"martin.aryee@gmail.com"
    }
}
task basecall_duplex  {
    input {
        String sample_id
        File fast5_archive
        String basecall_model
    }
    command <<<
        mkdir fast5s
        filetype=$(file ${fast5_archive})

        if [[ "$filetype" == *"gzip compressed data"* ]]; then
          echo "FAST5s appear to be compressed with gzip. Decompressing..."
          tar zxvf ${fast5_archive} -C fast5s
        fi

        if [[ "$filetype" == *"Zip archive data"* ]]; then
          echo "FAST5s appear to be compressed with zip. Decompressing..."
          unzip ${fast5_archive} -d fast5s
        fi
        
        # Convert FAST5 into POD5
        pod5 convert fast5 -r fast5s pod5s/reads.pod5 --threads 12 

        # Simplex call with --emit-moves
        dorado basecaller /dorado_models/${basecall_model} pod5s --emit-moves | samtools view -Sh > unmapped_reads_with_moves.bam

        # Identify potential pairs
        duplex_tools pair --output_dir ./pairs unmapped_reads_with_moves.bam
    
        # Stereo duplex basecall:
        dorado duplex /dorado_models/${basecall_model} pod5s --pairs pairs/pair_ids_filtered.txt | samtools view -Sh > ${sample_id}.duplex.bam
    >>>
    runtime {
    	gpuType: "nvidia-tesla-v100"
      	gpuCount: 1
      	nvidiaDriverVersion: "470.161.03"
      	zones: ["us-central1-a"] 
        docker: "us-central1-docker.pkg.dev/aryeelab/docker/dorado"
    }
    output {
        File duplex_bam = "${sample_id}.duplex.bam"
        File pairs = "pairs/pair_ids_filtered.txt"
    }
}
