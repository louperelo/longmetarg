FROM continuumio/miniconda3
LABEL author="Louisa Wessels Perelo" \
      description="Docker image containing all requirements for pipeline tools: blast 2.10, python=3.6"

COPY prokblast.yml /
RUN conda env create -f prokblast.yml && conda clean -a
RUN apt update && apt install -y procps
ENV PATH /opt/conda/envs/prokblast/bin:$PATH

# download BLAST prokaryont database 'ref_prok_rep_genomes'
RUN mkdir -p /blast_db
WORKDIR /blast_db
RUN update_blastdb.pl --decompress ref_prok_rep_genomes

#get BLAST taxid file 'taxdump.tar.gz'
#https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
RUN wget -O taxdump.tar.gz https://ftp.ncbi.nlm.nih.gov/pub/taxonomy 