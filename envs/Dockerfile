FROM debian:buster

RUN apt update && apt install -y \
    git \
    libxml2-dev \
    libxslt-dev \
    python3 \
    python3-biopython \
    python3-click \
    python3-jinja2 \
    python3-matplotlib \
    python3-pandas \
    python3-pdfkit \
    python3-pip \
    python3-pybedtools \
    python3-pysam \
    python3-pyvcf && \
    apt clean && apt purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*


RUN pip3 install svgutils crimson

RUN git clone https://github.com/brentp/align.git && \ 
    cd align && python3 setup.py install
