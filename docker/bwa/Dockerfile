#===============================================================================
#
#         FILE: Dockerfile
#    DEV USAGE: docker run -it -v /Users/Alec/Documents/Bioinformatics/MDV_Project/20190811_fastq2bam_steep/data:/Users/Alec/Documents/Bioinformatics/MDV_Project/20190811_fastq2bam_steep/data --name bwa ubuntu:latest
#        USAGE: docker image build -t steepale/bwa:1.0 . # local image build
#
#  DESCRIPTION:  This Dockerfile will build bwa
# REQUIREMENTS:  ---
#        NOTES:  bwa Github: https://github.com/lh3/bwa
#       AUTHOR:  Alec Steep, alec.steep@gmail.com
#  AFFILIATION:  Michigan State University (MSU), East Lansing, MI, United States
#				         USDA ARS Avian Disease and Oncology Lab (ADOL), East Lansing, MI, United States
#				         Technical University of Munich (TUM), Weihenstephan, Germany
#      VERSION:  1.0
#      CREATED:  2019.09.13
#     REVISION:  ---
#===============================================================================

# Pull the Ubuntu OS image
FROM ubuntu:18.04

# Change to base directory
WORKDIR /

# Install the bwa dependencies
# Note: use 'apt-cache policy <package-name>' to determine which package is installed
RUN apt-get update \
    && apt-get install -y sudo \
    && sudo apt-get install -y \
    wget \
    zip \
    cpanminus \
    build-essential \
    git

# Install zlib zip file in opt directory
WORKDIR /opt

# Install zlib
RUN wget http://zlib.net/zlib-1.2.11.tar.gz \
    && tar xvzf zlib-1.2.11.tar.gz \
    && rm zlib-1.2.11.tar.gz \
    && cd zlib-1.2.11 \
    && ./configure \
    && make test \
    && make install

# Install bwa zip file in opt directory
WORKDIR /opt
# Install bwa zip file in opt directory
RUN git clone https://github.com/lh3/bwa.git \
    && cd bwa \
    && make \
    && chmod 755 /opt/bwa/bwa \
    && ln -s /opt/bwa/bwa /usr/local/bin/bwa

# Back to the working directory
WORKDIR /

# This is the CMD command from the Ubuntu:18.04 image we FROM'ed
CMD ["/bin/bash"]