FROM ubuntu:22.04

RUN apt-get update && \
    apt-get -y install wget
RUN wget https://github.com/PacificBiosciences/pb-CpG-tools/releases/download/v2.3.2/pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu.tar.gz && \
    tar -xzf pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu.tar.gz
RUN ln -s /pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu/bin/aligned_bam_to_cpg_scores /usr/local/bin/ && \
    rm pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu.tar.gz