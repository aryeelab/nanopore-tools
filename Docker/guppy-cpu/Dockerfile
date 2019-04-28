FROM ubuntu:16.04

RUN apt-get update && apt-get install -y wget default-jdk apt-transport-https

RUN wget -O- https://mirror.oxfordnanoportal.com/apt/ont-repo.pub | apt-key add -  &&\
    echo "deb http://mirror.oxfordnanoportal.com/apt xenial-stable non-free" | tee /etc/apt/sources.list.d/nanoporetech.sources.list &&\
    apt-get update && \
    apt-get install -y ont-guppy-cpu=3.0.3-1~xenial


