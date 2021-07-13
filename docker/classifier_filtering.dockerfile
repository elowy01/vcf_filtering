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
				   libffi-dev \
				   libcurl4-openssl-dev \
				   libssl-dev \
                   tabix \
				   wget \
				   && apt-get clean


WORKDIR tmp/

#prepare Python
WORKDIR /tmp/python
RUN wget https://www.python.org/ftp/python/3.7.0/Python-3.7.0.tgz && tar xzvf Python-3.7.0.tgz
WORKDIR Python-3.7.0
RUN ./configure && make && make install
RUN ln -s /usr/local/bin/python3.7 /usr/bin/python && ln -s /usr/local/bin/pip3.7 /usr/bin/pip
RUN rm -rf /tmp/python/

# install BCFTools
RUN git clone --recurse-submodules git://github.com/samtools/htslib.git && git clone git://github.com/samtools/bcftools.git
WORKDIR bcftools
RUN make && make install
WORKDIR /tmp/
RUN rm -rf htslib && rm -rf bcftools

#install vt
WORKDIR /tmp/
RUN git clone https://github.com/atks/vt.git
WORKDIR vt/
RUN git submodule update --init --recursive 
RUN make
RUN cp vt /bin/
RUN rm -r /tmp/vt

#get the vcf_filtering repo
WORKDIR /lib/
RUN wget https://github.com/elowy01/vcf_filtering/archive/refs/tags/v1.0.2.tar.gz && tar -xvf v1.0.2.tar.gz
ENV PYTHONPATH=/lib/vcf_filtering-1.0.2/src/


#install Python libraries
RUN pip install pandas && pip install scikit-learn==0.20.3