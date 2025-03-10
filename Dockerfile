## Dockerfile
FROM ubuntu:22.04
MAINTAINER Amanda Cooksey
LABEL Description="AgBase Pathannotator"

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
    nano

RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py39_25.1.1-2-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh


# give write permissions to conda folder
RUN chmod 777 -R /opt/conda/

ENV PATH=$PATH:/opt/conda/bin

RUN conda config --add channels bioconda

RUN conda upgrade conda

RUN pip install pandas

# add hmmer and diamond

RUN conda install -c conda-forge -c bioconda agat

RUN conda install -c conda-forge -c bioconda orthofinder

RUN conda install --solver=classic -c conda-forge -c bioconda hmmer

ENV PATH /usr/bin/:$PATH

ADD pipeline/pathannotator.sh /usr/bin

ADD pipeline/pull_data.sh /usr/bin

ADD pipeline/merge_data.py /usr/bin

WORKDIR /usr/bin

RUN git clone https://github.com/takaram/kofam_scan.git

# Change the permissions and the path for the wrapper script
RUN chmod +x /usr/bin/pathannotator.sh

WORKDIR /root

RUN mkdir /workdir /data

RUN chmod a+w /workdir /data

# Entrypoint
ENTRYPOINT ["/usr/bin/pathannotator.sh"]

# Add path to working directory
WORKDIR /workdir
