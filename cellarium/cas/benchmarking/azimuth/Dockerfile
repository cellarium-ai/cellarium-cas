FROM satijalab/azimuth:0.5.0

RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 python3-pip \
    && rm -rf /var/lib/apt/lists/* \
    && ln -sf /usr/bin/python3 /usr/bin/python

RUN pip3 install anndata scipy

# Fix Matrix/SeuratObject version incompatibility: newer Matrix (>=1.6) removed the
# virtual class "mMatrix" that SeuratObject's LogMap depends on. Registering it in
# Rprofile.site ensures it exists before any package namespace is initialized.
RUN echo 'methods::setClass("mMatrix", contains = "VIRTUAL")' >> /usr/local/lib/R/etc/Rprofile.site

COPY helpers/convert_h5ad.py /scripts/convert_h5ad.py
COPY run_azimuth.R /scripts/run_azimuth.R
COPY run.sh /scripts/run.sh

ENTRYPOINT ["/scripts/run.sh"]
