FROM satijalab/azimuth:0.5.0

RUN apt-get update && apt-get install -y --no-install-recommends \
    wget ca-certificates bash \
    && rm -rf /var/lib/apt/lists/* \
    && wget -qO /tmp/miniforge.sh \
    https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh \
    && bash /tmp/miniforge.sh -b -p /opt/py311 \
    && rm /tmp/miniforge.sh \
    && /opt/py311/bin/conda install -y -c conda-forge \
    python=3.11 \
    "anndata>=0.11,<0.12" \
    scipy \
    pandas \
    h5py \
    && /opt/py311/bin/conda clean -afy \
    && echo 'methods::setClass("mMatrix", contains = "VIRTUAL")' >> /usr/local/lib/R/etc/Rprofile.site

COPY helpers/convert_h5ad.py /scripts/convert_h5ad.py
COPY run_azimuth.R /scripts/run_azimuth.R
COPY run.sh /scripts/run.sh

RUN chmod +x /scripts/run.sh \
    && useradd -m -u 1000 appuser \
    && chown -R appuser:appuser /scripts

USER appuser

ENTRYPOINT ["/scripts/run.sh"]