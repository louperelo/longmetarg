FROM continuumio/miniconda3
LABEL author="Louisa Wessels Perelo" \
      description="Docker image containing requirements for minimap2, samtools with CARD (rgi), python=3.6"

SHELL ["/bin/bash", "-c"]

COPY mapcard.yml /
#conda env with python 3.6, required by rgi (CARD database)
RUN conda env create -f mapcard.yml && conda clean -a
RUN apt update && apt install -y procps
ENV PATH /opt/conda/envs/mapcard/bin:$PATH
#RUN apt-get install -qq perl-doc

#from: https://github.com/arpcard/rgi/blob/master/Dockerfile
# download latest card database
RUN mkdir -p /card_db
WORKDIR /card_db
RUN wget -O data.tar.bz2 https://card.mcmaster.ca/latest/data && \
        mkdir -p canonical &&  tar xf data.tar.bz2 -C canonical && \
        rm data.tar.bz2

# install CARD database
RUN rgi card_annotation -i canonical/card.json > card_annotation.log 2>&1 && \
    rgi load -i canonical/card.json --card_annotation card_database_v*.fasta

#get BLAST taxid file 'taxdump.tar.gz'
#https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
RUN wget -O taxdump.tar.gz https://ftp.ncbi.nlm.nih.gov/pub/taxonomy 




