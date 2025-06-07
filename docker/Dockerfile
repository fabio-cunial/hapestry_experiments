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
        openjdk-17-jdk \
        openjdk-17-jre \
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
    
# BEDTOOLS
ARG bedtools_version=2.31.0
RUN wget https://github.com/arq5x/bedtools2/releases/download/v${bedtools_version}/bedtools.static \
    && mv bedtools.static bedtools \
    && chmod +x bedtools \
    && mv bedtools /usr/local/bin \
    && bedtools --help


# ----------------------------- NON-SV CALLERS ---------------------------------

ARG trgt_version=1.5.1
RUN wget https://github.com/PacificBiosciences/trgt/releases/download/v1.5.1/trgt-v${trgt_version}-x86_64-unknown-linux-gnu.tar.gz \
    && tar -zxvf trgt-v${trgt_version}-x86_64-unknown-linux-gnu.tar.gz \
    && mv trgt-v${trgt_version}-x86_64-unknown-linux-gnu/trgt . \
    && chmod +x trgt \
    && ./trgt --help


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


# ------------------------------- ASSEMBLERS -----------------------------------

ARG hifiasm_version=0.24.0
RUN wget https://github.com/chhylp123/hifiasm/archive/refs/tags/${hifiasm_version}.tar.gz \
    && tar -xzvf ${hifiasm_version}.tar.gz \
    && rm -f ${hifiasm_version}.tar.gz \
    && mv hifiasm-${hifiasm_version} hifiasm \
    && cd hifiasm \
    && make -j 16 \
    && ./hifiasm -h


# ------------------------------- SV CALLERS -----------------------------------

# SNIFFLES
ARG sniffles_version=2.5.3
RUN pip3 install sniffles==${sniffles_version} \
    && sniffles -h

# PBSV
ARG pbsv_version=2.10.0
RUN conda install -c bioconda pbsv=${pbsv_version} \
    && pbsv --help

# CUTESV
ARG cutesv_version=2.1.1
RUN wget https://github.com/tjiangHIT/cuteSV/archive/refs/tags/cuteSV-v${cutesv_version}.tar.gz \
    && tar xzf cuteSV-v${cutesv_version}.tar.gz \
    && rm -f cuteSV-v${cutesv_version}.tar.gz \
    && cd cuteSV-cuteSV-v${cutesv_version} \
    && python setup.py install \
    && cuteSV -h

# # SVISION
# ARG svision_version=2.4
# RUN wget https://github.com/songbowang125/SVision-pro/archive/refs/tags/v${svision_version}.tar.gz \
#     && tar xzf v${svision_version}.tar.gz \
#     && rm -f v${svision_version}.tar.gz \
#     && mv ./SVision-pro-2.4/ ./svision/ \
#     && cd ./svision/ \
#     && conda env create -f ./environment.yml \
#     && conda install mkl
# SHELL ["conda", "run", "-n", "svision-pro-env", "/bin/bash", "-c"]
# RUN cd ./svision/ \
#     && python setup.py install
# SHELL ["/bin/bash", "-c"]

# # SAWFISH
# ARG sawfish_version=0.12.9
# RUN wget https://github.com/PacificBiosciences/sawfish/releases/download/v0.12.9/sawfish-v${sawfish_version}-x86_64-unknown-linux-gnu.tar.gz \
#     && tar xzf sawfish-v${sawfish_version}-x86_64-unknown-linux-gnu.tar.gz \
#     && rm -f sawfish-v${sawfish_version}-x86_64-unknown-linux-gnu.tar.gz \
#     && mv ./sawfish-v${sawfish_version}-x86_64-unknown-linux-gnu/ ./sawfish/ \
#     && ./sawfish/bin/sawfish --help
    
# # NANOVAR
# ARG nanovar_version=1.8.3
# RUN wget https://github.com/benoukraflab/NanoVar/archive/refs/tags/v${nanovar_version}.tar.gz \
#     && tar xzf v${nanovar_version}.tar.gz \
#     && rm -f v${nanovar_version}.tar.gz \
#     && cd NanoVar-${nanovar_version} \
#     && pip install .
    
# SVIM
ARG svim_version=2.0.0
RUN wget https://github.com/eldariont/svim/archive/refs/tags/v${svim_version}.tar.gz \
    && tar xzf v${svim_version}.tar.gz \
    && rm -f v${svim_version}.tar.gz \
    && cd svim-${svim_version} \
    && pip install . \
    && svim --help
# Bypassing the final plotting, since it throws an error.
COPY ./svim /opt/conda/bin/svim

# # WTDBG2 (for debreak)
# ARG wtdbg2_version=2.5
# RUN wget https://github.com/ruanjue/wtdbg2/releases/download/v2.5/wtdbg-${wtdbg2_version}_x64_linux.tgz \
#     && tar xzf wtdbg-${wtdbg2_version}_x64_linux.tgz \
#     && rm -f wtdbg-${wtdbg2_version}_x64_linux.tgz \
#     && mv ./wtdbg-${wtdbg2_version}_x64_linux/ ./wtdbg2/ \
#     && ./wtdbg2/wtdbg2 --help

# # DEBREAK
# ARG debreak_version=1.2
# RUN wget https://github.com/Maggi-Chen/DeBreak/archive/refs/tags/v${debreak_version}.tar.gz \
#     && tar xzf v${debreak_version}.tar.gz \
#     && rm -f v${debreak_version}.tar.gz \
#     && mv ./DeBreak-${debreak_version}/ ./debreak/ \
#     && ./debreak/debreak --help

# # PAV (tried but failed to make the latest version work with the AoU WDLs)
# ARG pav_version=2.4.6
# ARG svpop_commit=968993798c770ce279c72a3004bbace1d2053b72
# RUN wget https://github.com/EichlerLab/pav/archive/refs/tags/v${pav_version}.tar.gz \
#     && tar xzf v${pav_version}.tar.gz \
#     && rm -f v${pav_version}.tar.gz \
#     && mv pav-${pav_version} pav \
#     && cd pav/dep \
#     && git clone --recursive https://github.com/EichlerLab/svpop.git \
#     && cd svpop \
#     && git checkout ${svpop_commit} \
#     && cd ../../.. \
#     && mkdir -p /opt/ \
#     && mv pav /opt/pav
# COPY ./pav_environment.yml .
# RUN conda env create -f pav_environment.yml \
#     && conda clean -a
# ENV PATH=/opt/conda/envs/lr-pav/bin/:${PATH}
# COPY ./pav_test_install.py /tmp/
# #SHELL ["conda", "run", "-n", "lr-pav", "/bin/bash", "-c"]
# #RUN python /tmp/pav_test_install.py


# ------------------------------- SV MERGERS -----------------------------------

# TRUVARI
ARG truvari_version=5.0.0
# ARG truvari_version=4.3.1   # OLD TRUVARI, JUST TO MEASURE CHANGES.
RUN pip3 install truvari==${truvari_version} \
    && truvari --help

# JASMINE
ARG jasmine_version=1.1.5
RUN conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && conda install jasminesv=${jasmine_version}

# PANPOP
ARG panpop_version=NC2024
RUN wget https://github.com/starskyzheng/panpop/archive/refs/tags/${panpop_version}.tar.gz \
    && tar xzf ${panpop_version}.tar.gz \
    && rm -f ${panpop_version}.tar.gz \
    && cd panpop-${panpop_version} \
    && conda env create -f conda.panpop.yaml
#SHELL ["conda", "run", "-n", "panpop", "/bin/bash", "-c"]

# # MAFFT
# RUN conda install conda-forge::mafft \
#     && mafft --help || echo "mafft returns an error code when called"

# # JASMINE  Some dependency fails now...
# COPY ./jasmine_environment.yml .
# RUN conda env create -f ./jasmine_environment.yml
# SHELL ["conda", "run", "-n", "jasmine", "/bin/bash", "-c"]
# Adding #!/bin/bash to the original shell script
# COPY ./jasmine /opt/conda/envs/jasmine/bin/jasmine
# SHELL ["chmod", "+x", "/opt/conda/envs/jasmine/bin/jasmine"]
# SHELL ["/bin/bash", "-c"]

# ------------------------------ SV GENOTYPERS ---------------------------------

ARG kanpig_version=1.0.2
RUN wget https://github.com/ACEnglish/kanpig/releases/download/v${kanpig_version}/kanpig-v${kanpig_version}-x86_64-unknown-linux-musl.tar.gz \
    && tar xzf kanpig-v${kanpig_version}-x86_64-unknown-linux-musl.tar.gz \
    && ls -laht kanpig-v${kanpig_version}-x86_64-unknown-linux-musl \
    && mv ./kanpig-v${kanpig_version}-x86_64-unknown-linux-musl/kanpig . \
    && chmod +x ./kanpig \
    && ./kanpig --help

# ---------------------------------- PHASERS -----------------------------------

ARG hiphase_version=1.4.5
RUN wget https://github.com/PacificBiosciences/HiPhase/releases/download/v${hiphase_version}/hiphase-v${hiphase_version}-x86_64-unknown-linux-gnu.tar.gz \
    && tar -xzf hiphase-v${hiphase_version}-x86_64-unknown-linux-gnu.tar.gz \
    && mv hiphase-v${hiphase_version}-x86_64-unknown-linux-gnu/hiphase ./hiphase \
    && chmod +x ./hiphase \
    && ./hiphase -h


# ----------------------------- CUSTOM SCRIPTS ---------------------------------
COPY *.java .
RUN javac *.java
COPY *.py .


# AWS
RUN wget https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip \
    && unzip awscli-exe-linux-x86_64.zip \
    && ./aws/install \
    && aws --version