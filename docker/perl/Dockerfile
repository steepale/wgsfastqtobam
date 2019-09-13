#===============================================================================
#
#         FILE: Dockerfile
#    DEV USAGE: docker run -it --name perl --rm -v /Users/Alec/Documents/Bioinformatics/MDV_Project/nf-core-wgsfastqtobam/data:/Users/Alec/Documents/Bioinformatics/MDV_Project/nf-core-wgsfastqtobam/data --name perl --rm ubuntu:18.04
#        USAGE: docker image build -t steepale/perl:5.26.1 . # local image build
#
#  DESCRIPTION:  This Dockerfile will create a perl installation on top of Ubuntu 18.04
# REQUIREMENTS:  ---
#        NOTES:  At time of construction, this image downloads perl version 5.26.1
#       AUTHOR:  Alec Steep, alec.steep@gmail.com
#  AFFILIATION:  Michigan State University (MSU), East Lansing, MI, United States
#				         USDA ARS Avian Disease and Oncology Lab (ADOL), East Lansing, MI, United States
#				         Technical University of Munich (TUM), Weihenstephan, Germany
#      VERSION:  5.26.1
#      CREATED:  2019.09.12
#     REVISION:  ---
#===============================================================================

# Base Ubuntu image
FROM ubuntu:18.04

# Set working dorectory
WORKDIR /

# Update aptitude and install some required packages
# a lot of them are required for Bio::DB::BigFile
RUN apt-get update \
    && apt-get -y install sudo \
    && sudo apt-get install -y \
    build-essential \
    cpanminus \
    perl \
    perl-base \
    wget \
    vim
# Note: use 'apt-cache policy <package-name>' to determine which package is installed

# Back to the working directory
WORKDIR /

# This is the CMD command from the Ubuntu:18.04 image we FROM'ed
CMD ["/bin/bash"]
