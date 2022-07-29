FROM ubuntu

RUN apt-get update && \
	DEBIAN_FRONTEND=noninteractive apt-get install -y \
		python3 \
		python-is-python3 \
		python3-numpy \
		python3-biopython \
		python3-seaborn \
		python3-sklearn \
		git

# Install nanopolish
RUN git clone https://github.com/al-mcintyre/mCaller.git