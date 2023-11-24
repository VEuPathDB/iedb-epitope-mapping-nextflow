FROM python:3.7
#FROM debian:stretch


ENV DEBIAN_FRONTEND=noninteractive
ENV biopython_version=1.81

USER root


RUN apt-get update 
RUN apt-get install -y git ant build-essential wget unzip bcftools tabix samtools perl \
default-jre unzip cpanminus bioperl libaio1 emacs libjson-perl libmodule-install-rdf-perl \
libdate-manip-perl libtext-csv-perl libstatistics-descriptive-perl libtree-dagnode-perl libxml-simple-perl && apt-get clean && apt-get purge && rm -rf /var/lib/apt/lists/* /tmp/*


RUN pip install --upgrade pip

WORKDIR /usr/local/

RUN pip install biopython==${biopython_version}


ADD /bin/* /usr/bin/
RUN cd /usr/bin \
  && chmod +x *.pl \
  && chmod +x *.py