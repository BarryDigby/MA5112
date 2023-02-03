FROM nfcore/base:1.14

COPY environment.yml /
RUN conda env create -f /environment.yml

ENV PATH /opt/conda/envs/week1/bin/:$PATH

