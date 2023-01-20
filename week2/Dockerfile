FROM nfcore/base:1.14

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

ENV PATH /opt/conda/envs/week2/bin/:$PATH
