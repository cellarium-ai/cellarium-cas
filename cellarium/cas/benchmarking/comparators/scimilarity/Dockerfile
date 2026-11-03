FROM python:3.12-slim

# g++ is required to build hnswlib from source (no wheel for this platform).
RUN apt-get update && apt-get install -y --no-install-recommends g++ \
    && rm -rf /var/lib/apt/lists/*

# CPU-only torch from the PyTorch CPU index (extra index: PyPI stays primary so
# torch's runtime deps resolve there). Installed before scimilarity so the
# resolver sees the torch requirement already satisfied and never pulls CUDA wheels.
RUN pip install --no-cache-dir --extra-index-url https://download.pytorch.org/whl/cpu \
    "torch==2.5.1+cpu" \
    && pip install --no-cache-dir \
    "scimilarity==0.4.1" \
    click

COPY run.py /scripts/run.py

RUN useradd -m -u 1000 appuser && chown -R appuser:appuser /scripts

USER appuser
ENTRYPOINT ["python", "/scripts/run.py"]
