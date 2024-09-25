FROM ubuntu:23.10

COPY requirement.r 

RUN apt-get update && apt install -y python2 sra-toolkit bcftools bwa wget unzip samtools \
tabix bzip2 openjdk-8-jdk git-lfs vim varscan python3  python-is-python3 bzip2 delly \
r-base wget curl xml2 gh openssl bedtools libxml2-dev libcurl4-openssl-dev libfontconfig1-dev \
libssl-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev

RUN git clone https://github.com/broadinstitute/gatk.git && cd gatk && git reset --hard 4.3.0.0 && ./gradlew bundle && ln -s /gatk/gatk /usr/bin/gatk

RUN wget https://github.com/Illumina/strelka/releases/download/v2.9.10/strelka-2.9.10.centos6_x86_64.tar.bz2 && tar -xvjf strelka-2.9.10.centos6_x86_64.tar.bz2 && \
rm strelka-2.9.10.centos6_x86_64.tar.bz2 && mv strelka-2.9.10.centos6_x86_64 strelka

RUN wget https://github.com/Illumina/manta/releases/download/v1.6.0/manta-1.6.0.centos6_x86_64.tar.bz2 && tar -xvjf manta-1.6.0.centos6_x86_64.tar.bz2 && \
rm manta-1.6.0.centos6_x86_64.tar.bz2 && mv manta-1.6.0.centos6_x86_64 manta

RUN wget https://cran.r-project.org/src/contrib/Archive/sequenza/sequenza_3.0.0.tar.gz
RUN Rscript requirement.r
RUN rm sequenza_3.0.0.tar.gz
