FROM rocker/rstudio:4.4.0

# Install required system packages
RUN apt-get update && apt-get install -y \
    openjdk-17-jdk \
    wget \
    unzip \
    python3 \
    python-is-python3 \
    libxml2-dev \
    libcurl4-openssl-dev \
    bzip2 \
    make \
    gcc \
    bwa \
    git \
    parallel \
    zlib1g-dev \
    libssl-dev \
    libffi-dev \
    libbz2-dev \
    liblzma-dev \
    libpcre3-dev \
    libicu-dev \
    tabix \
    perl \
    libncurses5-dev \
    libcurses-perl

# Install R packages in smaller groups
RUN R -e "install.packages(c('dplyr', 'httr', 'jsonlite'))"
RUN R -e "install.packages(c('optparse', 'purrr', 'stringr'))"
RUN R -e "install.packages(c('vcfR', 'BiocManager', 'remotes'))"
RUN R -e "install.packages(c('assertthat', 'xlsx', 'readr', 'knitr'))"
RUN R -e "remotes::install_version('matrixStats', version='1.4.1', repos='https://cloud.r-project.org')"

# Install Bioconductor packages in smaller groups
RUN R -e "BiocManager::install(version='3.20', ask=FALSE)"
RUN R -e "BiocManager::install('MatrixGenerics', ask=FALSE)"
RUN R -e "BiocManager::install(c('Biostrings', 'Rsamtools'), ask=FALSE)"
RUN R -e "BiocManager::install(c('VariantAnnotation', 'GenomicRanges'), ask=FALSE)"
RUN R -e "BiocManager::install('universalmotif', ask=FALSE)"
RUN R -e "BiocManager::install('regioneR', ask=FALSE)"

# Install Picard
RUN mkdir -p /opt/tools && \
    cd /opt/tools && \
    wget https://github.com/broadinstitute/picard/releases/download/3.1.1/picard.jar

# Install specific version of samtools
RUN cd /opt && \
    wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2 && \
    tar xjf samtools-1.10.tar.bz2 && \
    cd samtools-1.10 && \
    ./configure && \
    make && \
    make install && \
    ln -s /usr/local/bin/samtools /usr/bin/samtools  # Create symlink to expected location

# Install VarDict
RUN cd /opt && \
    git clone https://github.com/AstraZeneca-NGS/VarDictJava.git && \
    cd VarDictJava && \
    git clone https://github.com/AstraZeneca-NGS/VarDict && \
    cd VarDict && \
    chmod -R 755 . && \
    cp teststrandbias.R /usr/local/bin/ && \
    cp var2vcf_valid.pl /usr/local/bin/ && \
    chmod 755 /usr/local/bin/teststrandbias.R && \
    chmod 755 /usr/local/bin/var2vcf_valid.pl && \
    cd .. && \
    ./gradlew clean installDist && \
    chmod -R 755 build/install/VarDict && \
    ln -s /opt/VarDictJava/build/install/VarDict/bin/VarDict /usr/local/bin/ && \
    chmod 755 /usr/local/bin/VarDict

ENV VARDICT_HOME="/opt/VarDictJava/build/install/VarDict"
ENV PATH="$VARDICT_HOME/bin:/usr/local/bin:${PATH}"

# Install htslib for bgzip
RUN cd /opt && \
    wget https://github.com/samtools/htslib/releases/download/1.19/htslib-1.19.tar.bz2 && \
    tar xjf htslib-1.19.tar.bz2 && \
    cd htslib-1.19 && \
    ./configure && \
    make && \
    make install && \
    ldconfig

# Add PicardCommandLine wrapper script
RUN echo '#!/bin/sh\njava -jar /opt/tools/picard.jar SamToFastq "$@"' > /usr/local/bin/SamToFastq && \
    chmod +x /usr/local/bin/SamToFastq

# Set environment variables
ENV VARDICT_HOME="/opt/VarDictJava/build/install/VarDict"
ENV PATH="$VARDICT_HOME/bin:/usr/local/bin:${PATH}"
ENV PICARD_JAR="/opt/tools/picard.jar"

RUN chmod -R 755 /opt/VarDictJava/build/install/VarDict && \
    chmod -R 755 /usr/local/bin/VarDict && \
    chmod 755 /usr/local/bin/teststrandbias.R && \
    chmod 755 /usr/local/bin/var2vcf_valid.pl

RUN chown -R rstudio:rstudio /opt/VarDictJava

# Create reference directory and download/index reference genome
RUN mkdir -p /home/rstudio/reference/modified_reference

RUN mkdir -p /home/rstudio/results

# Download hg38 chromosome 7 reference
RUN cd /home/rstudio/reference && \
    wget -q https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr7.fa.gz && \
    gunzip chr7.fa.gz

# Index the reference genome with BWA and samtools
RUN cd /home/rstudio/reference && \
    bwa index chr7.fa && \
    samtools faidx chr7.fa

# Copy PMS2_vaR scripts
COPY PMS2_vaR/ /home/rstudio/PMS2_vaR/

# Create masked reference using R
RUN cd /home/rstudio/reference && \
    R -e "library(Biostrings); chr7_seq <- readDNAStringSet('chr7.fa'); pms2cl_start <- 6731863; pms2cl_end <- 6743066; chr7_string <- as.character(chr7_seq[[1]]); substr(chr7_string, pms2cl_start, pms2cl_end) <- paste(rep('N', pms2cl_end - pms2cl_start + 1), collapse=''); chr7_masked <- DNAStringSet(chr7_string); names(chr7_masked) <- 'chr7'; writeXStringSet(chr7_masked, filepath='modified_reference/chr7_masked.fa', width=60); cat('Masked reference created successfully\n')" && \
    cd modified_reference && \
    bwa index chr7_masked.fa && \
    samtools faidx chr7_masked.fa

# Set ownership for rstudio user
RUN chown -R rstudio:rstudio /home/rstudio/reference /home/rstudio/PMS2_vaR

WORKDIR /home/rstudio
