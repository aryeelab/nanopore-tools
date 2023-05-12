version 1.0

workflow demux {

    input {
        File infastq
    }

    call guppybarcoder {
        input:
            infastq=infastq
    }
    call makefastqs {
        input:
            fastq_gzs = guppybarcoder.fastq_gzs
    }
    output {
        File samplesheet=makefastqs.samplesheet
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
    guppy_barcoder -i ~{infastq} -s ./out --barcode_kits SQK-RBK114-24 --enable_trim_barcodes --compress_fastq
    >>>
    runtime {
		docker: "us-central1-docker.pkg.dev/aryeelab/docker/megalodon"
		memory: "64 GB"
		disks: "local-disk 1000 SSD"
		cpu: 12
    }
    output {
        Array[File] fastq_gzs = glob("./out/*/*.fq.gz")
    }
}
task makefastqs {
    input {
        Array[String] fastq_gzs
    }
    command <<<
        echo fastq_gz ~{sep=' ' fastq_gzs} >> samples_t.txt
        cat samples_t.txt | datamash --output-delimiter=',' -t ' ' transpose > samples.csv
    >>>
    runtime {
		docker: "us-central1-docker.pkg.dev/aryeelab/docker/megalodon"
		memory: "64 GB"
		disks: "local-disk 1000 SSD"
		cpu: 12
    }
    output {
        File samplesheet="samples.csv"
    }
}