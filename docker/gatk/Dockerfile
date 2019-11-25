#===============================================================================
#
#         FILE: Dockerfile
#    DEV USAGE: docker run -it -v /Users/Alec/Documents/Bioinformatics/MDV_Project/nf-core-wgsfastqtobam:/Users/Alec/Documents/Bioinformatics/MDV_Project/nf-core-wgsfastqtobam --name gatk steepale/ubuntu_jdk:1.8.0_212
#        USAGE: docker image build -t steepale/gatk: . # local image build
#
#  DESCRIPTION:  This Dockerfile will build gatk
# REQUIREMENTS:  ---
#        NOTES:  Picard Website: https://github.com/broadinstitute/picard/releases/tag/2.20.7
#       AUTHOR:  Alec Steep, alec.steep@gmail.com
#  AFFILIATION:  Michigan State University (MSU), East Lansing, MI, United States
#				         USDA ARS Avian Disease and Oncology Lab (ADOL), East Lansing, MI, United States
#				         Technical University of Munich (TUM), Weihenstephan, Germany
#      VERSION:  1.0
#      CREATED:  2019.09.13
#     REVISION:  ---
#===============================================================================

# Pull the Ubuntu OS image
FROM steepale/ubuntu_jdk:1.8.0_212

# Change to base directory
WORKDIR /

# Install the gatk dependencies
# Note: use 'apt-cache policy <package-name>' to determine which package is installed


# Install picard jar file in opt directory
WORKDIR /opt

# Install picard
RUN git clone https://github.com/broadinstitute/picard.git \
    && cd picard \
    && ./gradlew shadowJar \
    && chmod 755 /opt/picard/build/libs/picard.jar

RUN mkdir /opt/gatk3.5

# Install gatk zip file in opt directory
WORKDIR /opt/gatk3.5

# Install gatk zip file in opt directory
# GATK3.5 needed to be interactively installed: https://software.broadinstitute.org/gatk/download/archive
RUN bzip2 -d GenomeAnalysisTK-3.5-0-g36282e4.tar.bz2 \
    && tar xvf GenomeAnalysisTK-3.5-0-g36282e4.tar \
    && rm GenomeAnalysisTK-3.5-0-g36282e4.tar \
    && chmod 755 /opt/gatk3.5/GenomeAnalysisTK.jar

# Back to the working directory
WORKDIR /

# This is the CMD command from the Ubuntu:18.04 image we FROM'ed
CMD ["/bin/bash"]