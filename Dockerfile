FROM ubuntu:18.04

LABEL description="methyl-data analysis pipeline"

RUN apt-get update && apt-get install -y \
      git \
      g++ \
      build-essential \
      openjdk-8-jdk \
      openjdk-8-jre \
      perl \
      wget \
      curl \
      cpanminus \
      zlib1g \
      zlib1g-dev \
      python-pip \
      python3-pip \
      python \
			vim \
			tar

RUN apt-get update -y && \
				apt-get upgrade -y && \
				apt-get dist-upgrade -y && \
				apt-get install build-essential software-properties-common -y && \
				apt-get install libcurl4-openssl-dev -y && \
				apt-get install libssl-dev -y && \
				apt-get install libbz2-dev -y && \
				apt-get update -y

RUN cpanm \
      Math::Gauss \
      Math::Round \
      Bio::Trace::ABIF \
      Spreadsheet::Write \
			Statistics::Basic \
      Parallel::ForkManager \
	  Sort::Key

ENV CONDA_DIR /opt/conda
ENV PATH $CONDA_DIR/bin:$PATH

<<<<<<< HEAD
=======




#RUN wget --quiet --no-check-certificate https://repo.anaconda.com/archive/Anaconda3-2019.03-Linux-x86_64.sh && \
#    echo "45c851b7497cc14d5ca060064394569f724b67d9b5f98a926ed49b834a6bb73a	*Anaconda3-2019.03-Linux-x86_64.sh" | sha256sum -c - && \
#    /bin/bash /Anaconda3-2019.03-Linux-x86_64.sh -f -b -p $CONDA_DIR && \
#    rm Anaconda3-2019.03-Linux-x86_64.sh && \
#    echo export PATH=$CONDA_DIR/bin:'$PATH' > /etc/profile.d/conda.sh


>>>>>>> 8bd99478acddf4ab6f72ca7b810f24915461cd8a
RUN wget -q --no-check-certificate https://repo.anaconda.com/archive/Anaconda3-2022.05-Linux-x86_64.sh
RUN echo "a7c0afe862f6ea19a596801fc138bde0463abcbce1b753e8d5c474b506a2db2d *Anaconda3-2022.05-Linux-x86_64.sh" | sha256sum -c - && \
    /bin/bash /Anaconda3-2022.05-Linux-x86_64.sh -f -b -p $CONDA_DIR && \
    rm Anaconda3-2022.05-Linux-x86_64.sh && \
    echo export PATH=$CONDA_DIR/bin:'$PATH' > /etc/profile.d/conda.sh

ENV TZ=Europe/Moscow
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
RUN apt remove r-base
RUN add-apt-repository -r 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'
RUN apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/'


## R 4.1.3
RUN apt-get update -y
RUN apt-get install gdebi-core -y
RUN export R_VERSION=4.1.3
RUN curl -O https://cdn.rstudio.com/r/ubuntu-1804/pkgs/r-4.1.3_1_amd64.deb
RUN gdebi -nq /r-4.1.3_1_amd64.deb
RUN ln -s /opt/R/4.1.3/bin/R /usr/local/bin/R
RUN ln -s /opt/R/4.1.3/bin/Rscript /usr/local/bin/Rscript


RUN apt-get install -y pigz
RUN apt-get install libtbb2

RUN conda install -c bioconda bedtools 
RUN conda install -c bioconda ucsc-fasize
RUN conda install -c bioconda ucsc-wigtobigwig
RUN conda install -c bioconda fastqc
RUN conda install -c bioconda samtools

<<<<<<< HEAD
## bs2
RUN pip2 install --upgrade setuptools
RUN pip2 install pysam
RUN conda install -c bioconda bowtie


## bismark
#RUN conda install -c bioconda bowtie2
RUN git clone https://github.com/BenLangmead/bowtie2.git
WORKDIR bowtie2
RUN make
WORKDIR cd /
ENV PATH /bowtie2:$PATH

#RUN conda install -c bioconda bismark
#RUN conda install -c bioconda/label/cf201901 bismark

=======
RUN conda install -c bioconda bowtie

###
RUN conda install -c bioconda/label/broken bowtie2

RUN conda install -c bioconda bismark
>>>>>>> 8bd99478acddf4ab6f72ca7b810f24915461cd8a
RUN pip install multiqc
RUN pip install cutadapt

## trimgalore
RUN curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.6.tar.gz -o trim_galore.tar.gz
RUN tar xvzf trim_galore.tar.gz
ENV PATH /TrimGalore-0.6.6:$PATH

<<<<<<< HEAD
=======

RUN git clone "https://github.com/BSSeeker/BSseeker2.git"
ENV BS2 /BSseeker2
ENV PATH $BS2:$PATH
>>>>>>> 8bd99478acddf4ab6f72ca7b810f24915461cd8a

RUN mkdir -p /msPIPE
COPY lib /msPIPE/lib
COPY bin /msPIPE/bin
COPY msPIPE.py /msPIPE/
RUN Rscript /msPIPE/bin/script/Package_install.R


RUN wget https://github.com/bedops/bedops/releases/download/v2.4.39/bedops_linux_x86_64-v2.4.39.tar.bz2
RUN tar jxvf bedops_linux_x86_64-v2.4.39.tar.bz2
ENV PATH /msPIPE:$PATH
ENV PATH $CONDA_DIR/bin/cutadapt:$PATH


RUN mkdir -p /work_dir
WORKDIR /work_dir

