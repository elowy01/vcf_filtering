#parent image
FROM ubuntu:latest

LABEL maintainer="ernestolowy@gmail.com"
LABEL description="Dockerfile used to build an image used in genomic variant filtering.The filtering used is based on a supervised logistic regression classifier" 

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get -qq install git \
    				   build-essential \
				   autoconf \
				   zlib1g-dev \
				   libbz2-dev \
				   liblzma-dev \
				   libhts-dev  \
                   libvcflib-tools \
				   libvcflib-dev \
				   python3 \
				   python3-pip \
				   libcurl4-openssl-dev \
				   libssl-dev \
                   tabix \
				   && apt-get clean
				   
WORKDIR tmp/

#prepare Python
RUN ln -s /usr/bin/python3 /usr/bin/python

# install BCFTools
RUN git clone --recurse-submodules git://github.com/samtools/htslib.git && git clone git://github.com/samtools/bcftools.git
WORKDIR bcftools
RUN make && make install
WORKDIR /tmp/
RUN rm -rf htslib && rm -rf bcftools

#install igsr-analysis libraries
WORKDIR /lib
RUN git clone https://github.com/igsr/igsr_analysis.git
ENV PYTHONPATH=/lib/igsr_analysis
ENV PATH=/bin/:/lib/igsr_analysis/scripts/VCF/QC/:${PATH}

#install vt
WORKDIR /tmp/
RUN git clone https://github.com/atks/vt.git
WORKDIR vt/
RUN git submodule update --init --recursive 
RUN make
RUN cp vt /bin/
RUN rm -r /tmp/vt
WORKDIR /root/

RUN pip install pandas sklearn