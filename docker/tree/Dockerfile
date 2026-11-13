FROM mambaorg/micromamba:2.0.5@sha256:8d0d4cad074afe3ac93b58dce88bd2659153d3fa6de3a6968f501710e8d3c832
RUN micromamba create -y -n tree -c conda-forge -c bioconda iqtree=3.1.4 clonalframeml=1.20 mafft=7.526 \
      pal2nal=14.1 python=3.12 numpy=2.1.3 && \
    micromamba clean -a -y
COPY --chown=$MAMBA_USER:$MAMBA_USER pyproject.toml /tmp/plock/
COPY --chown=$MAMBA_USER:$MAMBA_USER src /tmp/plock/src
RUN micromamba run -n tree pip install --no-cache-dir --no-deps /tmp/plock && rm -rf /tmp/plock
RUN micromamba run -n tree iqtree3 --version && micromamba run -n tree ClonalFrameML -version && \
    micromamba run -n tree mafft --version && micromamba run -n tree sh -c "command -v pal2nal.pl" && \
    micromamba run -n tree plock_helpers --help >/dev/null
