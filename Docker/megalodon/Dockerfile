FROM nvidia/cuda:11.7.0-base-ubuntu20.04

RUN apt-get update && \
  		apt-get -y install python3-pip wget libbz2-dev liblzma-dev git

# Set up time zone for non-interactive tzdata install
RUN export DEBIAN_FRONTEND=noninteractive && \
  		ln -fs /usr/share/zoneinfo/America/New_York /etc/localtime

# Install guppy
RUN wget -nv https://cdn.oxfordnanoportal.com/software/analysis/ont_guppy_5.0.16-1~focal_amd64.deb && \
  		apt -y install ./ont_guppy_5.0.16-1~focal_amd64.deb 
  
# Install megalodon and pod5 tools
# Note that numpy<1.24 is needed since megalodon uses 'np.int' rather than 'int' (deprecated/removed in 1.24)
RUN pip install numpy==1.23 && \
    pip install megalodon pod5

RUN pip install ont-pyguppy-client-lib==5.0.16
