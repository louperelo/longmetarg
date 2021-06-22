FROM continuumio/miniconda3
LABEL author="Louisa Wessels Perelo" \
      description="Docker image containing requirements for PlasFlow"

#WORKDIR $PWD plasflow
#VOLUME $HOME:$HOME

COPY plasflow.yml /
RUN conda env create -f plasflow.yml && conda clean -a
RUN apt update && apt install -y procps
ENV PATH /opt/conda/envs/plasflow/bin:$PATH

#copy folder models into plasflow/bin
RUN mkdir /opt/conda/envs/plasflow/bin/models
COPY models /opt/conda/envs/plasflow/bin/models/


