FROM continuumio/miniconda3
LABEL author="Louisa Wessels Perelo" \
      description="Docker image containing requirements for pipeline tools: MetaFlye, python=3.8"

COPY flye.yml /
RUN conda env create -f flye.yml && conda clean -a
RUN apt update && apt install -y procps
ENV PATH /opt/conda/envs/flye/bin:$PATH



