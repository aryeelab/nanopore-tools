# standalone docker image 
FROM ubuntu:18.04

RUN apt-get update && apt-get -y upgrade
RUN apt-get install -y libc6-dev build-essential tcsh devscripts debhelper wget tar git
RUN apt-get clean

RUN git clone https://github.com/bedops/bedops.git
# RUN wget https://github.com/bedops/bedops/releases/download/v2.4.41/bedops_linux_x86_64-v2.4.41.tar.bz2
# RUN tar jxvf bedops_linux_x86_64-v2.4.41.tar.bz2
# copy the source context into the local image
#  note: make sure .dockerignore is up to date
COPY . /bedops

# build and install bedops into system path
WORKDIR /bedops
RUN make -j `nproc` && make install BINDIR=/usr/bin
WORKDIR /
RUN rm -rf /bedops