# Base Image
FROM python:3.7.3

# Metadata
LABEL base.image="translocationfinder:latest"
LABEL version="1"
LABEL software="TranslocationFinder"
LABEL software.version="latest"
LABEL description=""
LABEL tags="Translocation Sequencing NGS"

# Maintainer
MAINTAINER Dave Lab <lab.dave@gmail.com>

# update the OS related packages and other dependencies
RUN apt-get update -y &&\
    apt-get install samtools bedtools libz-dev

# install required python dependencies
RUN pip install numpy scipy pandas pybedtools pysam

# make directory to store tools such as QCParser
RUN mkdir tools

# get the QCParser from GitHub
RUN git clone --branch master https://github.com/labdave/TranslocationFinder.git /tools/TranslocationFinder

ENV PATH /tools/TranslocationFinder:$PATH

CMD ["python", ""]
