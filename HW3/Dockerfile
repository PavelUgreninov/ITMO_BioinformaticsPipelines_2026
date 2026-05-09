FROM debian:bullseye-slim

RUN apt-get update --fix-missing && \
    apt-get install -y wget bzip2 ca-certificates libglib2.0-0 libxext6 libsm6 libxrender1 git && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN wget --quiet https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O /tmp/miniforge.sh && \
    /bin/bash /tmp/miniforge.sh -b -p /opt/conda && \
    rm /tmp/miniforge.sh

ENV PATH="/opt/conda/bin:${PATH}"

COPY env.yaml /tmp/env.yaml

RUN mamba env update -n base -f /tmp/env.yaml && \
    mamba install -y -c bioconda bcftools htslib && \
    mamba clean -afy

WORKDIR /data

CMD ["/bin/bash"]