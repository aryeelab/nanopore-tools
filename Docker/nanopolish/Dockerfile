# From the nanopolish Dockerfile (https://github.com/jts/nanopolish/blob/master/Dockerfile)
# We modify it by pinning the version number and installing JDK for jar 
FROM centos:7
WORKDIR /
RUN yum group install "Development Tools" -y
RUN yum install git wget tar zlib-devel -y
RUN yum install java-1.7.0-openjdk-devel -y
RUN git clone --recursive https://github.com/jts/nanopolish.git
WORKDIR /nanopolish
RUN git checkout v0.11.1 && make all
CMD ./nanopolish