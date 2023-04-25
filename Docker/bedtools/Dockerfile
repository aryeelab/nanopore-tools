FROM ubuntu:18.04

RUN apt-get update && apt-get -y upgrade &&\
    apt-get install -y wget libcurl4-openssl-dev libkrb5-3 bedtools &&\
    apt-get clean && apt-get purge &&\
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
WORKDIR /usr/bin
RUN wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig && \
  chmod +x bedGraphToBigWig