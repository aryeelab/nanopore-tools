FROM nvidia/cuda:9.0-cudnn7-devel-ubuntu16.04

RUN apt-get update && apt-get install -y libcurl4-openssl-dev \
                                         libssl-dev \
                                         libhdf5-cpp-11 \
                                         libzmq5 \
                                         libboost-atomic1.58.0 \
                                         libboost-chrono1.58.0 \
                                         libboost-date-time1.58.0 \
                                         libboost-filesystem1.58.0 \
                                         libboost-program-options1.58.0 \
                                         libboost-regex1.58.0 \
                                         libboost-system1.58.0 \
                                         libboost-log1.58.0 \
                                         wget \
                                         default-jdk

RUN cd /tmp &&\
    wget -q https://mirror.oxfordnanoportal.com/software/analysis/ont_guppy_3.0.3-1~xenial_amd64.deb &&\
    dpkg -i --ignore-depends=nvidia-384,libcuda1-384 /tmp/ont_guppy_3.0.3-1~xenial_amd64.deb


