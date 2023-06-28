version 1.0

workflow demux {

    input {
        File infastq
        File reads
        File pythonscript
    }

    call guppybarcoder {
        input:
            infastq=infastq
    }
    call makefast5s {
        input:
            reads=reads,
            pythonscript=pythonscript,
            barcoding_summary = guppybarcoder.barcoding_summary
    }
    call makesheet {
        input:
            fastqs = guppybarcoder.fastqs,
            fast5_tars = makefast5s.fast5_tars
    }
    output {
        File samplesheet=makesheet.samplesheet
        Array[File] fast5_tars = makefast5s.fast5_tars
        Array[File] fastqs = guppybarcoder.fastqs
    }

    meta {
        author: "Mark Soliman"
        email:"msoliman@ds.dfci.harvard.edu"
    }
}
task guppybarcoder  {
    input {
        File infastq
    }
    command <<<
    mkdir out
    mkdir in
    cp ~{infastq} ./in
    guppy_barcoder -i ./in -s ./out --barcode_kits SQK-RBK114-24 #--enable_trim_barcodes
    for f in ./out/*; do
        dir_name="${f##*/}"
        cat "$f"/*.fastq > ./out/"$dir_name".fastq
        # for file in "$f"/*.fastq; do
        #     mv "$file" "$f/$dir_name-${file##*/}"
        # done
    done
    column=$(head -n 2 ./out/barcoding_summary.txt | tail -n 1 | awk -F$'\t' 'BEGIN{search="unclassified|barcode"} { for (i=1; i<=NF; i++) { if ($i ~ search) print i } }')
    cat ./out/barcoding_summary.txt | cut -f 1,${column} >> ./out/twocolumnsummary.tsv
    >>>
    runtime {
		docker: "us-central1-docker.pkg.dev/aryeelab/docker/megalodon"
		memory: "64 GB"
		disks: "local-disk 1000 SSD"
		cpu: 12
    }
    output {
        File barcoding_summary = "./out/twocolumnsummary.tsv"
        Array[File] fastqs = glob("./out/*.fastq")
    }
}
task makefast5s {
    input {
        File reads
        File pythonscript
        File barcoding_summary
    }
    command <<<
        filetype=$(file ~{reads})
        mkdir in
        if [[ ~{reads} == *.pod5 ]]; then
            pod5 convert to_fast5 ~{reads} --output ./in
            else
            tar xvzf ~{reads} -C ./in
        fi
        mkdir out
        column_name=$(awk -F'\t' 'NR==1 {print $2}' ~{barcoding_summary})
        python3 ~{pythonscript} --input ./in --save_path ./out --summary_file ~{barcoding_summary} --demultiplex_column ${column_name}
        for f in ./out/*; do tar czvf "$f.tar.gz" "$f"/*.fast5; done
    >>>
    runtime {
		docker: "us-central1-docker.pkg.dev/aryeelab/docker/ontfast5api"
		memory: "64 GB"
		disks: "local-disk 1000 SSD"
		cpu: 12
    }
    output {
        Array[File] fast5_tars= glob("./out/*.tar.gz")
    }
}
task makesheet {
    input {
        Array[String] fastqs
        Array[String] fast5_tars
    }
    command <<<
        echo fastq_gz ~{sep=' ' fastqs} >> samples_t.txt
        echo fast5_tar ~{sep=' ' fast5_tars} >> samples_t.txt
        cat samples_t.txt | datamash --output-delimiter=',' -t ' ' transpose > samples.csv
    >>>
    runtime {
        docker: "quay.io/aryeelab/nanopore-util"
		memory: "64 GB"
		disks: "local-disk 1000 SSD"
		cpu: 12
    }
    output {
        File samplesheet= "samples.csv"
    }
}