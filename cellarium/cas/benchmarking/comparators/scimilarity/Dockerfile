FROM ghcr.io/genentech/scimilarity@sha256:ca30629a3737e318230400e8d4d3d26615a0d7370d926562cabce622db99d248

RUN python3 -m pip install --no-cache-dir "scimilarity==0.4.1" \
    && ln -sf "$(command -v python3)" /usr/local/bin/python

COPY run.py /scripts/run.py

RUN useradd -m -u 1000 appuser && chown -R appuser:appuser /scripts

USER appuser
ENTRYPOINT ["python", "/scripts/run.py"]
