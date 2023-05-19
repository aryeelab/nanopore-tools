version 1.0

workflow demux {

    input {
        File infastq
        File reads
    }

    call guppybarcoder {
        input:
            infastq=infastq
    }
    call makefast5s {
        input:
            reads=reads,
            barcoding_summary = guppybarcoder.barcoding_summary
    }
    call makesheet {
        input:
            fastq_gzs = guppybarcoder.fastq_gzs,
            fast5_tars = makefast5s.fast5_tars
    }
    output {
        File samplesheet=makesheet.samplesheet
        Array[File] fast5_tars = makefast5s.fast5_tars
        Array[File] fastq_gzs = guppybarcoder.fastq_gzs
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
    guppy_barcoder -i ./in -s ./out --barcode_kits SQK-RBK114-24 --enable_trim_barcodes --compress_fastq
    >>>
    runtime {
		docker: "us-central1-docker.pkg.dev/aryeelab/docker/megalodon"
		memory: "64 GB"
		disks: "local-disk 1000 SSD"
		cpu: 12
    }
    output {
        File barcoding_summary = "./out/barcoding_summary.txt"
        Array[File] fastq_gzs = glob("./out/*/*.fastq.gz")
    }
}
task makefast5s {
    input {
        File reads
        File barcoding_summary
    }
    command <<<
        mkdir in
        tar xvzf ~{reads} -C ./in
        mkdir out
        demux_fast5 --input ./in --save_path ./out --summary_file ~{barcoding_summary}
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
        Array[String] fastq_gzs
        Array[String] fast5_tars
    }
    command <<<
        echo fastq_gz ~{sep=' ' fastq_gzs} >> samples_t.txt
        echo fast5_tar ~{sep=' ' fast5_tars} >> samples_t.txt
        cat samples_t.txt | datamash --output-delimiter=',' -t ' ' transpose > samples.csv
    >>>
    runtime {
        docker: "us-central1-docker.pkg.dev/aryeelab/docker/ontfast5api"
		memory: "64 GB"
		disks: "local-disk 1000 SSD"
		cpu: 12
    }
    output {
        File samplesheet= "samples.csv"
    }
}