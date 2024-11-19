FROM ubuntu:23.10

COPY requirement.r /

RUN apt-get update && apt install -y sra-toolkit bcftools bwa wget unzip samtools \
tabix bzip2 openjdk-8-jdk git-lfs vim varscan python3  python-is-python3 bzip2 delly \
r-base wget curl xml2 gh openssl bedtools libxml2-dev libcurl4-openssl-dev libfontconfig1-dev make \
libssl-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev autoconf libboost-all-dev cmake


RUN git clone --recursive https://github.com/dellytools/sansa.git && cd sansa && make all && ln -s /sansa/src/sansa /usr/bin/sansa

RUN wget https://www.python.org/ftp/python/2.7/Python-2.7.tgz && \
tar -xvzf Python-2.7.tgz && cd Python-2.7 && ./configure && make && ln -s /Python-2.7/python /usr/bin/python2.7 && ln -s /Python-2.7/python /usr/bin/python2

RUN git clone https://github.com/broadinstitute/gatk.git && cd gatk && git reset --hard 4.3.0.0 && ./gradlew bundle && ln -s /gatk/gatk /usr/bin/gatk

RUN wget https://github.com/Illumina/strelka/releases/download/v2.9.10/strelka-2.9.10.centos6_x86_64.tar.bz2 && tar -xvjf strelka-2.9.10.centos6_x86_64.tar.bz2 && \
rm strelka-2.9.10.centos6_x86_64.tar.bz2 && mv strelka-2.9.10.centos6_x86_64 strelka

RUN wget https://github.com/Illumina/manta/releases/download/v1.6.0/manta-1.6.0.centos6_x86_64.tar.bz2 && tar -xvjf manta-1.6.0.centos6_x86_64.tar.bz2 && \
rm manta-1.6.0.centos6_x86_64.tar.bz2 && mv manta-1.6.0.centos6_x86_64 manta

RUN wget https://cran.r-project.org/src/contrib/Archive/sequenza/sequenza_3.0.0.tar.gz
RUN Rscript requirement.r
RUN rm sequenza_3.0.0.tar.gz
