#===============================================================================
#
#         FILE: Dockerfile
#    DEV USAGE: docker run -it -v /Users/Alec/Documents/Bioinformatics/MDV_Project/nf-core-wgsfastqtobam:/Users/Alec/Documents/Bioinformatics/MDV_Project/nf-core-wgsfastqtobam --name ubuntu_jdk --rm ubuntu_jdk:11.0.4
#        USAGE: docker image build -t steepale/ubuntu_jdk:11.0.4 . # local image build
#
#  DESCRIPTION:  This Dockerfile will build ubuntu_jdk
# REQUIREMENTS:  ---
#        NOTES:  
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

# Install the java development kit dependencies
# Note: use 'apt-cache policy <package-name>' to determine which package is installed
RUN apt-get update \
    && apt-get install -y sudo \
    && sudo apt-get install -y \
    wget \
    zip \
    cpanminus \
    build-essential \
    git \
    default-jdk

# Back to the working directory
WORKDIR /

# This is the CMD command from the Ubuntu:18.04 image we FROM'ed
CMD ["/bin/bash"]

