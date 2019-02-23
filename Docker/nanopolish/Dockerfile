# From the nanopolish Dockerfile (https://github.com/jts/nanopolish/blob/master/Dockerfile)
# We only modify it by pinning the version number
FROM centos:7
WORKDIR /
RUN yum group install "Development Tools" -y
RUN yum install git wget tar zlib-devel -y
RUN git clone --recursive https://github.com/jts/nanopolish.git
WORKDIR /nanopolish
RUN git checkout v0.11.0 && make all
CMD ./nanopolish