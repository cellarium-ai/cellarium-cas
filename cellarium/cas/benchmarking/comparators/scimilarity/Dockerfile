FROM python:3.12-slim

# g++ builds hnswlib from source (PyPI ships an sdist only).
RUN apt-get update && apt-get install -y --no-install-recommends g++ \
    && rm -rf /var/lib/apt/lists/*

# Only the packages we install directly are pinned; scimilarity's transitive
# deps resolve freely. torch is in the same command so the resolver keeps the
# CPU build (+cpu local tag exists only on the CPU index; PyPI stays primary).
# --only-binary :all: everywhere except hnswlib, which has no wheel.
RUN pip install --no-cache-dir --only-binary :all: --no-binary hnswlib \
    --extra-index-url https://download.pytorch.org/whl/cpu \
    "torch==2.5.1+cpu" \
    "scimilarity==0.4.1" \
    "click==8.5.0"

COPY run.py /scripts/run.py

RUN useradd -m -u 1000 appuser && chown -R appuser:appuser /scripts

USER appuser
ENTRYPOINT ["python", "/scripts/run.py"]
