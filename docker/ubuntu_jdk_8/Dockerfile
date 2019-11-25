#===============================================================================
#
#         FILE: Dockerfile
#    DEV USAGE: docker run -it -v /Users/Alec/Documents/Bioinformatics/MDV_Project/nf-core-wgsfastqtobam:/Users/Alec/Documents/Bioinformatics/MDV_Project/nf-core-wgsfastqtobam --name ubuntu_jdk steepale/ubuntu_jdk:1.8.0_212
#        USAGE: docker image build -t steepale/ubuntu_jdk:1.8.0_212 . # local image build
#
#  DESCRIPTION:  This Dockerfile will build steepale/ubuntu_jdk:1.8.0_212
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
FROM steepale/ubuntu_base:18.04

# Change to base directory
WORKDIR /

# Install the java development kit dependencies
# Note: use 'apt-cache policy <package-name>' to determine which package is installed
RUN sudo apt-get update \
    && sudo apt-get install -y \
    software-properties-common \
    vim

# Install in opt directory
WORKDIR /opt

# The following file needed to be interactively downloaded from oracle. The following command will not work.
RUN wget https://www.oracle.com/technetwork/java/javase/downloads/java-archive-javase8u211-later-5573849.html \
    && tar xvzf jdk-8u212-linux-x64.tar.gz \
    && rm jdk-8u212-linux-x64.tar.gz \
    && sudo ln -s /opt/jdk1.8.0_212/bin/java /usr/bin/java \
    && echo 'JAVA_HOME="/usr/bin/java"' >> /etc/environment \
    source /etc/environment \
    chmod 755 /usr/bin/java

# Back to the working directory
WORKDIR /

# This is the CMD command from the Ubuntu:18.04 image we FROM'ed
CMD ["/bin/bash"]

