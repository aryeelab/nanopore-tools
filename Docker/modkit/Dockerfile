FROM centos:7
RUN yum -y install wget && mkdir /modkit
RUN wget -O - https://github.com/nanoporetech/modkit/releases/download/v0.2.5-rc2/modkit_v0.2.5-rc2_centos7_x86_64.tar.gz | tar -xzvf - -C /modkit
ENV PATH="${PATH}:/modkit/dist/"