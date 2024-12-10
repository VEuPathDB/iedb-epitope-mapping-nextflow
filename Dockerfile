FROM python:3.7
#FROM debian:stretch

ENV DEBIAN_FRONTEND=noninteractive
ENV biopython_version=1.81

RUN apt-get update 
RUN apt-get install -y git ant build-essential wget perl \
default-jre unzip cpanminus bioperl libaio1 libjson-perl libmodule-install-rdf-perl \
libdate-manip-perl libtext-csv-perl libstatistics-descriptive-perl libtree-dagnode-perl libxml-simple-perl && apt-get clean && apt-get purge && rm -rf /var/lib/apt/lists/* /tmp/*

RUN pip install --upgrade pip

RUN pip install biopython==${biopython_version}

RUN pip install pepmatch==1.0.5

WORKDIR /data
