FROM debian:buster

RUN apt update && apt install -y \
    git \
    libxml2-dev \
    libxslt-dev \
    libqt5core5a \
    python3 \
    python3-biopython \
    python3-jinja2 \
    python3-matplotlib \
    python3-pandas \
    python3-pdfkit \
    python3-pip \
    python3-pybedtools \
    python3-pysam \
    python3-pyvcf \
    r-bioc-gviz \
    r-cran-getopt && \
    apt clean && apt purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN apt update && apt install -y wget

RUN wget https://github.com/wkhtmltopdf/wkhtmltopdf/releases/download/0.12.5/wkhtmltox_0.12.5-1.buster_amd64.deb \
    && apt install -y ./wkhtmltox_0.12.5-1.buster_amd64.deb

RUN pip3 install svgutils crimson

RUN git clone https://github.com/brentp/align.git && \ 
    cd align && python3 setup.py install

COPY Dockerfile .
