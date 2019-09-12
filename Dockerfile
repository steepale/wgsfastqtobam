FROM nfcore/base
LABEL authors="Alec Steep & Hongen Xu" \
      description="Docker image containing all requirements for nf-core/wgsfastqtobam pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-wgsfastqtobam-1.0dev/bin:$PATH
