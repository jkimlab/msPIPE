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

RUN wget --quiet --no-check-certificate https://repo.anaconda.com/archive/Anaconda3-2019.03-Linux-x86_64.sh && \
    echo "45c851b7497cc14d5ca060064394569f724b67d9b5f98a926ed49b834a6bb73a *Anaconda3-2019.03-Linux-x86_64.sh" | sha256sum -c - && \
    /bin/bash /Anaconda3-2019.03-Linux-x86_64.sh -f -b -p $CONDA_DIR && \
    rm Anaconda3-2019.03-Linux-x86_64.sh && \
    echo export PATH=$CONDA_DIR/bin:'$PATH' > /etc/profile.d/conda.sh

ENV TZ=Europe/Moscow
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
RUN apt remove r-base
RUN add-apt-repository -r 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/'
RUN apt-get update -y
RUN apt-get dist-upgrade -y
RUN apt-get install r-base r-base-core r-recommended -y
RUN apt-get install -y pigz
RUN apt-get install libtbb2

RUN conda install -c bioconda bedtools
RUN conda install -c bioconda ucsc-fasize
RUN conda install -c bioconda ucsc-wigtobigwig
RUN conda install -c bioconda cutadapt
RUN conda install -c bioconda fastqc
RUN conda install -c bioconda samtools
RUN conda install -c bioconda bowtie2

RUN conda install -c bioconda trim-galore
RUN conda install -c bioconda bismark
#RUN conda install -c bioconda multiqc
RUN pip install multiqc

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

COPY params* /work_dir


