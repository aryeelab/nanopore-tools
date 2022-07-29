FROM ubuntu
	
RUN apt-get update && apt-get install -y \
		build-essential \
    	git \
    	samtools \
    	zlib1g-dev \
    	&& rm -rf /var/lib/apt/lists/*
	
RUN git clone https://github.com/lh3/minimap2

RUN cd minimap2 && \
	if [ "$(uname -m)" = "aarch64" ]; then \
    	make arm_neon=1 aarch64=1; \
    else \
        make; \
    fi
