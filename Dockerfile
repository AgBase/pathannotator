## Dockerfile
FROM ubuntu:22.04
MAINTAINER Amanda Cooksey
LABEL Description="AgBase Pathannotator"

ENV DEBIAN_FRONTEND=noninteractive

# Install all the updates and download dependencies
RUN apt-get update && \
    apt-get install -y \
    git \
    wget \
    gzip \
    parallel \
    python3 \
    ruby \
    tar \
    nano \
    python3-scipy \
    python3-sklearn \
    python3-numpy \
    python3-biopython \
    liblist-moreutils-perl \
    libtry-tiny-perl \
    libbio-perl-perl \
    libclone-perl \
    libgraph-perl \
    liblwp-useragent-determined-perl \
    libstatistics-r-perl \
    libcarp-clan-perl \
    libsort-naturally-perl \
    libfile-share-perl \
    libfile-sharedir-perl \
    libfile-sharedir-install-perl \
    libyaml-perl \
    liblwp-protocol-https-perl \
    libterm-progressbar-perl

RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py39_25.1.1-2-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh


# give write permissions to conda folder
RUN chmod 777 -R /opt/conda/

ENV PATH=$PATH:/opt/conda/bin

RUN conda config --add channels bioconda

RUN conda upgrade conda

#RUN conda install -c conda-forge -c bioconda cd-hit

RUN conda install -c conda-forge -c bioconda agat

RUN conda install -c conda-forge -c bioconda orthofinder=3.0.1b1-0

RUN conda install --solver=classic -c conda-forge -c bioconda hmmer

RUN conda install --solver=classic -c conda-forge -c bioconda pandas

RUN conda install -c conda-forge -c bioconda seqkit

ENV PERL5LIB=$PERL5LIB:/opt/conda/pkgs:/opt/conda/pkgs/agat-1.4.2-pl5321hdfd78af_1/lib/perl5/site_perl/:/opt/conda/bin

ENV PATH /usr/bin/:$PATH

ADD pipeline/pathannotator.sh /usr/bin

ADD pipeline/pull_data.sh /usr/bin

ADD pipeline/merge_data.py /usr/bin

ADD pipeline/build_ref_set.sh /usr/bin

RUN mkdir /AGAT /OF

ADD pipeline/agat_config.yaml /AGAT

COPY pipeline/ref_set.tgz /OF

WORKDIR /usr/bin

RUN git clone https://github.com/takaram/kofam_scan.git

RUN wget http://github.com/bbuchfink/diamond/releases/download/v2.1.11/diamond-linux64.tar.gz && tar xzf diamond-linux64.tar.gz

# Change the permissions and the path for the wrapper script
RUN chmod +x /usr/bin/pathannotator.sh build_ref_set.sh

WORKDIR /root

RUN mkdir /workdir /data 

RUN chmod a+w /workdir /data /OF /AGAT

# Entrypoint
ENTRYPOINT ["/usr/bin/pathannotator.sh"]

# Add path to working directory
WORKDIR /workdir

