FROM ubuntu

RUN apt-get update && \
	apt-get install -y \
		build-essential \
		cmake \
		git \
		libbz2-dev \
		libhdf5-dev \
		libzstd-dev \
		seqkit \
		tar \
		wget \
		zlib1g-dev


# Install nanopolish
COPY config.guess.patch /tmp/
RUN git clone --recursive https://github.com/jts/nanopolish.git && \
	cd nanopolish && \
	git apply /tmp/config.guess.patch && \
	make all

# Install VBZ compression support
RUN cd / && \
	git clone https://github.com/nanoporetech/vbz_compression.git && \
	cd vbz_compression && \
	git submodule update --init && \
	mkdir build && \
	cd build && \
	cmake -D CMAKE_BUILD_TYPE=Release -D ENABLE_CONAN=OFF -D ENABLE_PERF_TESTING=OFF -D ENABLE_PYTHON=OFF .. && \
	make -j6 && \
	make install

ENV HDF5_PLUGIN_PATH /usr/local/hdf5/lib/plugin

ENV PATH="/nanopolish:${PATH}"
ARG VERSION=__specify_in_build_arg__
RUN echo $VERSION > /VERSION

WORKDIR /nanopolish

ENTRYPOINT [ "nanopolish" ]