FROM mambaorg/micromamba:2.0.5@sha256:8d0d4cad074afe3ac93b58dce88bd2659153d3fa6de3a6968f501710e8d3c832
RUN micromamba create -y -n panaroo -c conda-forge -c bioconda panaroo=1.8.1 && \
    micromamba clean -a -y
RUN micromamba run -n panaroo panaroo --help >/dev/null
