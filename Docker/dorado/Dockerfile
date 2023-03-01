FROM nvidia/cuda:12.0.1-devel-ubuntu22.04

# Install basic packages
RUN apt-get update && \
    apt-get -y install  file \
                        libbz2-dev \
                        liblzma-dev \
                        python3-pip \
                        samtools \
                        unzip \
                        wget

# Install pod5 tools and duplex tools
RUN pip install pod5 duplex_tools

# Install Dorado
RUN wget -q https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.2.1-linux-x64.tar.gz && \
    tar zxf dorado-0.2.1-linux-x64.tar.gz && \
    ln -s /dorado-0.2.1-linux-x64/bin/dorado /usr/local/bin/ && \
    rm dorado-0.2.1-linux-x64.tar.gz

# Install basecalling models
RUN mkdir /dorado_models && \
    cd /dorado_models && \
    dorado download 
