FROM continuumio/miniconda3:24.1.2-0
MAINTAINER Fabio Cunial
ARG work_dir=/hapestry
WORKDIR ${work_dir}


# -------------------------------- BASICS --------------------------------------

# OS
RUN apt-get -qqy update --fix-missing \
    && apt-get -qqy dist-upgrade \
    && apt-get install -y --no-install-recommends \
        zlib1g-dev \
        liblzma-dev \
        libbz2-dev \
        libdeflate-dev \
        libreadline-dev \
        libsqlite3-dev \
        libssl-dev \
        libcurl4-openssl-dev \
        libncurses5-dev \
        libncursesw5-dev \
        libffi-dev \
        liblzma-dev \
        libopenblas-dev \
        apt-transport-https \
        gawk \
        ca-certificates \
        tree \
        gnupg \
        ssh \
        time \
        curl \
        wget \
        autotools-dev \
        autoconf \
        automake \
        make \
        cmake \
        gcc \
        g++ \
        gfortran \
        build-essential \
        git \
        bc \
        python3-pip \
        xz-utils \
        tk-dev \
        python-dev \
        bsdmainutils \
        default-jre \
        default-jdk \
        unzip \
        python3-pycurl \
        python3-edlib \
        pigz

# GSUTIL
ARG gcloud_version=405.0.1
RUN pip3 uninstall -y crcmod && pip3 install --no-cache-dir -U crcmod
RUN wget https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-cli-${gcloud_version}-linux-x86_64.tar.gz \
    && tar -xf google-cloud-cli-${gcloud_version}-linux-x86_64.tar.gz \
    && rm -f google-cloud-cli-${gcloud_version}-linux-x86_64.tar.gz \
    && yes | ./google-cloud-sdk/install.sh
ENV PATH=${work_dir}/google-cloud-sdk/bin:${PATH}

# HTSLIB
ARG htslib_version=1.19.1
RUN wget https://github.com/samtools/htslib/releases/download/${htslib_version}/htslib-${htslib_version}.tar.bz2 \
    && tar xjf htslib-${htslib_version}.tar.bz2 \
    && rm htslib-${htslib_version}.tar.bz2 \
    && cd htslib-${htslib_version} \
    && ./configure \
    && make -j 16 \
    && make install \
    && cd ${work_dir} \
    && rm -rf htslib-${htslib_version} \
    && bgzip --help

# SAMTOOLS
ARG samtools_version=1.19.2
RUN wget https://github.com/samtools/samtools/releases/download/${samtools_version}/samtools-${samtools_version}.tar.bz2 \
    && tar xjf samtools-${samtools_version}.tar.bz2 \
    && rm samtools-${samtools_version}.tar.bz2 \
    && cd samtools-${samtools_version} \
    && ./configure --without-curses \
    && make -j 16 \
    && make install \
    && cd ${work_dir} \
    && rm -rf samtools-${samtools_version} \
    && samtools --help

# BCFTOOLS
ARG bcftools_version=1.21
RUN wget https://github.com/samtools/bcftools/releases/download/${bcftools_version}/bcftools-${bcftools_version}.tar.bz2 \
    && tar xjf bcftools-${bcftools_version}.tar.bz2 \
    && rm bcftools-${bcftools_version}.tar.bz2 \
    && cd bcftools-${bcftools_version} \
    && ./configure --without-curses \
    && make -j 16 \
    && make install \
    && cd ${work_dir} \
    && bcftools --help
    
# SEQKIT
ARG seqkit_version=2.8.2
RUN wget https://github.com/shenwei356/seqkit/releases/download/v2.8.2/seqkit_linux_amd64.tar.gz \
    && tar -xzvf seqkit_linux_amd64.tar.gz \
    && rm -f seqkit_linux_amd64.tar.gz \
    && chmod +x seqkit \
    && ./seqkit --help
    
# # BEDTOOLS
# ARG bedtools_version=2.31.0
# RUN wget https://github.com/arq5x/bedtools2/releases/download/v${bedtools_version}/bedtools.static \
#     && mv bedtools.static bedtools \
#     && chmod +x bedtools \
#     && mv bedtools /usr/local/bin \
#     && bedtools --help


# ------------------------------- ALIGNMENT ------------------------------------

# MINIMAP2
ARG minimap2_version=2.28
RUN wget https://github.com/lh3/minimap2/releases/download/v${minimap2_version}/minimap2-${minimap2_version}_x64-linux.tar.bz2 \
    && tar xjf minimap2-${minimap2_version}_x64-linux.tar.bz2 \
    && rm -f minimap2-${minimap2_version}_x64-linux.tar.bz2 \
    && mv ./minimap2-${minimap2_version}_x64-linux/ ./minimap2/ \
    && ./minimap2/minimap2 -h && echo 0 || echo 1

# PBMM2
ARG pbmm2_version=1.16.0
RUN conda install -c bioconda pbmm2=${pbmm2_version} \
    && pbmm2 --help


# ----------------------------- VARIANT CALLING --------------------------------

# # SNIFFLES
# ARG sniffles_version=2.2
# RUN pip3 install sniffles==${sniffles_version} \
#     && sniffles -h
#
# # PBSV
# ARG pbsv_version=2.9.0
# RUN conda install -c bioconda pbsv=${pbsv_version} \
#     && pbsv --help

# # CUTESV
# RUN pip install setuptools==58 scikit-learn \
#     && git clone https://github.com/tjiangHIT/cuteSV.git \
#     && cd cuteSV \
#     && python setup.py install \
#     && cuteSV -h



# ------------------------------- SV MERGING -----------------------------------

# # TRUVARI
# ARG truvari_version=4.2.2
# RUN pip3 install truvari==${truvari_version} \
#     && truvari --help

# # MAFFT
# RUN conda install conda-forge::mafft \
#     && mafft --help || echo "mafft returns an error code when called"

# # JASMINE  Some dependency fails now...
# COPY ./jasmine_environment.yml .
# RUN conda env create -f ./jasmine_environment.yml
# SHELL ["conda", "run", "-n", "jasmine", "/bin/bash", "-c"]
# # Adding #!/bin/bash to the original shell script
# COPY ./jasmine /opt/conda/envs/jasmine/bin/jasmine
# SHELL ["chmod", "+x", "/opt/conda/envs/jasmine/bin/jasmine"]
# SHELL ["/bin/bash", "-c"]


# ----------------------------- CUSTOM SCRIPTS ---------------------------------
COPY *.java .
RUN javac *.java
