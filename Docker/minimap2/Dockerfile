FROM ubuntu:20.04
	
RUN apt-get update && apt-get install -y build-essential wget libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev libcurl3-dev \
		bzip2 \
		curl &&\
		apt-get autoclean \
    	&& rm -rf /var/lib/apt/lists/*
	
RUN curl -L https://github.com/lh3/minimap2/releases/download/v2.26/minimap2-2.26_x64-linux.tar.bz2 | tar -jxvf - && \
 mkdir /data

RUN wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2 && \
	tar jxf samtools-1.16.1.tar.bz2 -C /minimap2-2.26_x64-linux && \
	rm samtools-1.16.1.tar.bz2 && \
	cd minimap2-2.26_x64-linux && \
	cp -a samtools-1.16.1/. ./ && \
	rm -r samtools-1.16.1/ && \
	./configure --prefix $(pwd) && \
	make
ENV PATH="${PATH}:/minimap2-2.26_x64-linux"
WORKDIR /data
