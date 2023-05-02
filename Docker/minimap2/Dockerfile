FROM ubuntu:20.04
	
RUN apt-get update && apt-get install -y \
    	samtools \
    	zlib1g-dev \
		bzip2 \
		curl &&\
		apt-get autoclean \
    	&& rm -rf /var/lib/apt/lists/*
	
RUN curl -L https://github.com/lh3/minimap2/releases/download/v2.26/minimap2-2.26_x64-linux.tar.bz2 | tar -jxvf - && \
 mkdir /data

ENV PATH="${PATH}:/minimap2-2.26_x64-linux"
WORKDIR /data
