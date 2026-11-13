FROM --platform=linux/arm64 mambaorg/micromamba:2.0.5@sha256:8d0d4cad074afe3ac93b58dce88bd2659153d3fa6de3a6968f501710e8d3c832 AS diamond-arm64
RUN micromamba create -y -n diamond -c conda-forge -c bioconda diamond=2.2.8

FROM --platform=linux/amd64 mambaorg/micromamba:2.0.5@sha256:8d0d4cad074afe3ac93b58dce88bd2659153d3fa6de3a6968f501710e8d3c832
RUN micromamba create -y -n ggcaller -c conda-forge -c bioconda ggcaller=1.5.2
RUN micromamba create -y -n bakta -c conda-forge -c bioconda bakta=1.12.2 diamond=2.2.8
RUN micromamba clean -a -y
COPY --from=diamond-arm64 /opt/conda/envs/diamond /opt/diamond-arm64/env
COPY --from=diamond-arm64 /lib/aarch64-linux-gnu /opt/diamond-arm64/system-lib
RUN mv /opt/conda/envs/bakta/bin/diamond /opt/conda/envs/bakta/bin/diamond-x86_64
COPY --chmod=755 docker/annotation/diamond /opt/conda/envs/bakta/bin/diamond
COPY --chown=$MAMBA_USER:$MAMBA_USER pyproject.toml /tmp/plock/
COPY --chown=$MAMBA_USER:$MAMBA_USER src /tmp/plock/src
RUN /opt/conda/envs/bakta/bin/pip install --no-cache-dir --no-deps /tmp/plock && rm -rf /tmp/plock
USER root
RUN ln -s /opt/conda/envs/bakta/bin/plock_helpers /usr/local/bin/plock_helpers
USER $MAMBA_USER
RUN micromamba run -n ggcaller ggcaller --help >/dev/null && \
    micromamba run -n bakta bakta_proteins --help >/dev/null && \
    plock_helpers --help >/dev/null
