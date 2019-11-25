#===============================================================================
#
#         FILE: Dockerfile
#    DEV USAGE: docker run -it -v /Users/Alec/Documents/Bioinformatics/MDV_Project/20190811_fastq2bam_steep/data:/Users/Alec/Documents/Bioinformatics/MDV_Project/20190811_fastq2bam_steep/data --name trimmomatic --rm steepale/trimmomatic/1.0
#        USAGE: docker image build -t steepale/trimmomatic:1.0 . # local image build
#
#  DESCRIPTION:  This Dockerfile will build trimmomatic
# REQUIREMENTS:  ---
#        NOTES:  Trimmomatic Website: http://www.usadellab.org/cms/?page=trimmomatic
#       AUTHOR:  Alec Steep, alec.steep@gmail.com
#  AFFILIATION:  Michigan State University (MSU), East Lansing, MI, United States
#				         USDA ARS Avian Disease and Oncology Lab (ADOL), East Lansing, MI, United States
#				         Technical University of Munich (TUM), Weihenstephan, Germany
#      VERSION:  1.0
#      CREATED:  2019.09.12
#     REVISION:  ---
#===============================================================================

# Pull the Ubuntu OS image
FROM ubuntu:18.04

# Change to base directory
WORKDIR /

# Install the Trimmomatic dependencies
# Note: use 'apt-cache policy <package-name>' to determine which package is installed
RUN apt-get update \
    && apt-get install -y sudo \
    && sudo apt-get install -y \
    default-jre \
    wget \
    vim \
    less \
    zip \
    cpanminus

# Install trimmomatic file in opt directory
WORKDIR /opt

# Install Trimmomatic, adjust permissions, & create sym link in local bin
RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip \
    && unzip Trimmomatic-0.39.zip \
    && rm Trimmomatic-0.39.zip \
    && chmod 755 /opt/Trimmomatic-0.39/trimmomatic-0.39.jar \
    && ln -s /opt/Trimmomatic-0.39/trimmomatic-0.39.jar /usr/local/bin/trimmomatic.jar

# Back to the working directory
WORKDIR /

# This is the CMD command from the Ubuntu:18.04 image we FROM'ed
CMD ["/bin/bash"]