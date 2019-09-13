#===============================================================================
#
#         FILE: Dockerfile
#    DEV USAGE: docker run -it -v /Users/Alec/Documents/Bioinformatics/MDV_Project/nf-core-wgsfastqtobam/data:/Users/Alec/Documents/Bioinformatics/MDV_Project/nf-core-wgsfastqtobam/data --name fastqc --rm steepale/fastqc:1.0
#        USAGE: docker image build -t steepale/fastqc:1.0 . # local image build
#
#  DESCRIPTION:  This Dockerfile will build FastQC
# REQUIREMENTS:  ---
#        NOTES:  FastQC Website: www.bioinformatics.babraham.ac.uk/projects/fastqc/
#       AUTHOR:  Alec Steep, alec.steep@gmail.com
#  AFFILIATION:  Michigan State University (MSU), East Lansing, MI, United States
#				         USDA ARS Avian Disease and Oncology Lab (ADOL), East Lansing, MI, United States
#				         Technical University of Munich (TUM), Weihenstephan, Germany
#      VERSION:  1.0
#      CREATED:  2019.08.11
#     REVISION:  ---
#===============================================================================

# Pull the Ubuntu OS image
FROM ubuntu:18.04

# Change to base directory
WORKDIR /

# Install the FastQC dependencies
# Note: use 'apt-cache policy <package-name>' to determine which package is installed
RUN apt-get update \
    && apt-get install -y sudo \
    && sudo apt-get install -y default-jdk \
    wget \
    zip \
    cpanminus

# Configure java
RUN echo 'JAVA_HOME="/usr/bin/java"' >> /etc/environment \
    && source /etc/environment \
    && chmod 755 /usr/bin/java

# Install FastQC zip file in opt directory
WORKDIR /opt

# Install FastQC, adjust permissions, & create sym link in local bin
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip \
    && unzip fastqc_v0.11.8.zip \
    && rm fastqc_v0.11.8.zip \
    && chmod 755 /opt/FastQC/fastqc \
    && mkdir /root/bin \
    && sudo ln -s /opt/FastQC/fastqc /usr/local/bin/fastqc 

# Back to the working directory
WORKDIR /

# This is the CMD command from the Ubuntu:18.04 image we FROM'ed
CMD ["/bin/bash"]