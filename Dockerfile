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
#agat    liblist-moreutils-perl \
#agat    libbio-perl-perl \
#agat    libclone-perl \
#agat    libgraph-perl \ 
#agat    liblwp-useragent-determined-perl \
#agat    libstatistics-r-perl \
#agat    libcarp-clan-perl \
#agat    libsort-naturally-perl \
#agat    libfile-share-perl \
#agat    libfile-sharedir-install-perl \
#agat    libyaml-perl \
#agat    liblwp-protocol-https-perl \
#agat    libfile-sharedir-perl \
#agat    libmoose-perl \
#agat    libterm-progressbar-perl \
#agat    libdevel-cover-perl \
    python3-scipy \
    python3-sklearn \
    python3-numpy \
    python3-biopython

RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py39_25.1.1-2-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh


# give write permissions to conda folder
RUN chmod 777 -R /opt/conda/

ENV PATH=$PATH:/opt/conda/bin

RUN conda config --add channels bioconda

RUN conda upgrade conda

#RUN pip install pandas 

# add hmmer and AGAT and orthofinder

#RUN conda install -c conda-forge -c bioconda agat

#RUN conda install -c conda-forge -c bioconda gffread

RUN conda install -c conda-forge -c bioconda cd-hit

RUN conda install -c conda-forge -c bioconda orthofinder

RUN conda install --solver=classic -c conda-forge -c bioconda hmmer

RUN conda install --solver=classic -c conda-forge -c bioconda pandas

ENV PERL5LIB=$PERL5LIB:/opt/conda/pkgs:/opt/conda/pkgs/agat-1.4.2-pl5321hdfd78af_1/lib/perl5/site_perl/:/opt/conda/bin

ENV PATH /usr/bin/:$PATH

ADD pipeline/pathannotator.sh /usr/bin

ADD pipeline/pull_data.sh /usr/bin

ADD pipeline/merge_data.py /usr/bin

#ADD pipeline/agat_config.yaml /usr/bin

WORKDIR /usr/bin

RUN git clone https://github.com/takaram/kofam_scan.git

RUN wget http://github.com/bbuchfink/diamond/releases/download/v2.1.11/diamond-linux64.tar.gz && tar xzf diamond-linux64.tar.gz

# Change the permissions and the path for the wrapper script
RUN chmod +x /usr/bin/pathannotator.sh

WORKDIR /root

RUN mkdir /workdir /data

RUN chmod a+w /workdir /data

# Entrypoint
ENTRYPOINT ["/usr/bin/pathannotator.sh"]

# Add path to working directory
WORKDIR /workdir

