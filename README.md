# Oxford nanopore tools

## Installing dependencies

- Install Homebrew

https://brew.sh

- Install cromwell

```
brew install cromwell
```

- Install Docker

https://hub.docker.com/editions/community/docker-ce-desktop-mac


## mcaller m6A pipeline

See https://github.com/al-mcintyre/mcaller


    #SAMPLE_NAME="20220107-ctcf-hia5"
    WORK_DIR="${PWD}/work/20220107_CTCF_hia5_pilot/20220701_1749_MC-111988_0_FAR26597_5eec18a4"    

### Download the genome

    GENOME_URL="ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz"
    curl ${GENOME_URL} -o ${WORK_DIR}/genome.fna.gz

### Build the Nanopolish Docker image

    cd ..
    git clone --recursive https://github.com/aryeelab/nanopolish.git
    cd nanopolish
    docker build -t nanopolish .
    cd ../nanopore-mods

### Build the minimap2 Docker image

    cd Docker/minimap2
    docker build -t minimap2 .
    cd ../..

### Build the mcaller Docker image

    cd Docker/mcaller
    docker build -t mcaller .
    cd ../..

    


### Index reads and raw signal data

    docker run -v ${WORK_DIR}/:/work/ --rm -it nanopolish /bin/bash
        cd work
        #zcat fastq_pass/*.fastq.gz > reads.fastq
        zcat fastq_pass/*_0.fastq.gz > reads.fastq
        seqkit stats reads.fastq > reads.fastq.stats
        time /nanopolish/nanopolish index -d fast5_pass -s sequencing_summary*.txt reads.fastq 2> stderr.log   
        cat stderr.log  | grep 'num reads' | sed 's/\[readdb\] //' | tr , '\n' | sed 's/^ //' > nanopolish-index.stats.txt

### Align reads to the genome. 

    docker run -v ${WORK_DIR}/:/work/ --rm -it minimap2 /bin/bash
        cd /work
        time /minimap2/minimap2 -a genome.fna.gz reads.fastq | \
            samtools view -Sb - | \
            samtools sort -o reads.sorted.bam 
        samtools index reads.sorted.bam 
        samtools flagstat reads.sorted.bam > reads.sorted.bam.flagstat 


### Align events to the genome

    docker run -v ${WORK_DIR}/:/work/ --rm -it nanopolish /bin/bash
        cd /work
        gunzip genome.fna.gz
        /nanopolish/nanopolish eventalign -t 6 --scale-events -n -r reads.fastq -b reads.sorted.bam -g genome.fna 1> reads.eventalign.tsv 2> stderr.log
        cat stderr.log | sed 's/\[post-run summary\] //' | tr , '\n' | sed 's/^ //' > eventalign.stats.txt

### Detect m6A

    docker run -v ${WORK_DIR}/:/work/ --rm -it mcaller /bin/bash
        cd /mCaller
        gunzip genome.fna.gz        
        # 5 hours CTCF
        time ./mCaller.py -m GATC -t 5 -r /work/genome.fna -d r95_twobase_model_NN_6_m6A.pkl -e /work/reads.eventalign.tsv -f /work/reads.fastq -b A  
        time ./make_bed.py -f /work/reads.eventalign.diffs.6 -d 1 -t 0


### Make m6A bigwigs
    
    docker run -v ${WORK_DIR}/:/work/ --rm -it 4dndcic/4dn-bedgraphtobigwig /bin/bash
        
        cd /work
        CHROM_SIZES="hg38.chrom.sizes"
        wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/${CHROM_SIZES}
        
        # A few reads are misassigned to the wrong alt contig? Temp fix:
        bedClip -verbose=2 reads.methylation.summary.bed ${CHROM_SIZES} reads.methylation.summary.clipped.bed        

        cat reads.methylation.summary.clipped.bed | cut -f1,2,3,5 | sort -k1,1 -k2,2n > m6a.ratio.bedgraph
        cat reads.methylation.summary.clipped.bed | cut -f1,2,3,7 | sort -k1,1 -k2,2n > m6a.coverage.bedgraph
        
        bedGraphToBigWig  m6a.ratio.bedgraph ${CHROM_SIZES} m6a.ratio.bw
        bedGraphToBigWig  m6a.coverage.bedgraph ${CHROM_SIZES} m6a.coverage.bw




## Basecalling/mC pipeline

### Quick Start

```
cromwell run -i test_data/test-run-1.cpu.json preprocess_flowcell.wdl 
```

You can include the path to a custom Cromwell config file by something like:

```
JAVA_OPTS="-Dconfig.file=/Users/maryee/cromwell/cromwell.conf"
cromwell run -i test_data/test-run-1.json preprocess_flowcell.wdl 
```
(See https://cromwell.readthedocs.io/en/latest/tutorials/ConfigurationFiles/)

This can be used to specify a cache database, for example:

```
docker run -p 3306:3306 --name cromwell_mysql -e MYSQL_ROOT_PASSWORD=pass123 -e MYSQL_DATABASE=cromwell -e MYSQL_USER=cromwell -e MYSQL_PASSWORD=pass123 -d mysql/mysql-server:5.5

```
(See https://cromwell.readthedocs.io/en/latest/tutorials/PersistentServer/)

### Testing within a docker image

e.g.
```
docker run --rm -it aryeelab/guppy
```

### Run wdl in erisone with singularity
java -jar -Dconfig.file=lsf.bash.conf /data/molpath/software/cromwell-39.jar run preprocess_flowcell_singularity.wdl -i test_data/test-run-1.json

## Visualizing the workflow graph

You can use `womtool` (part of cromwell) to output a workflow graph in `.dot` format:

```
womtool graph preprocess_flowcell.wdl > preprocess_flowcell.dot
```

This graph can be edited if necessary (such as to label edges with inputs/outputs), and then visualized with graphviz (`brew install graphviz`):

```
dot preprocess_flowcell.dot -Tpng -o preprocess_flowcell.png
```

The workflow graph below is produced in this way.


## Preprocess Flowcell workflow

![alt text](preprocess_flowcell.png "preprocess_flowcell.wdl DAG")

