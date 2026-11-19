FROM python:3.12-slim

# CPU-only; celltypist's declared dependencies (numpy, pandas, scikit-learn, scanpy, ...) resolve upstream.
RUN pip install --no-cache-dir --only-binary :all: \
    "click==8.5.0" \
    "celltypist==1.7.1"

COPY run.py /scripts/run.py

RUN useradd -m -u 1000 appuser && chown -R appuser:appuser /scripts

# Models are mounted read-only and passed via --model-path; never downloaded. CellTypist creates
# its cache directories under $HOME on import.
ENV HOME=/home/appuser

USER appuser
ENTRYPOINT ["python", "/scripts/run.py"]
