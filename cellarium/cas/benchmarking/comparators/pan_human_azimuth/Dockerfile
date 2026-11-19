FROM python:3.11-slim

# Only the packages we install directly are pinned; panhumanpy's transitive deps
# (TensorFlow 2.17, scikit-learn 1.6.0, ...) resolve upstream.
RUN pip install --no-cache-dir --only-binary :all: \
    "click==8.5.0" \
    "panhumanpy==1.0.0"

COPY run.py /scripts/run.py

RUN useradd -m -u 1000 appuser && chown -R appuser:appuser /scripts \
    && install -d -o appuser -g appuser /home/appuser/.cache/panhumanpy

# CPU-only. Model weights download to $HOME/.cache/panhumanpy; mount it to persist.
ENV HOME=/home/appuser \
    CUDA_VISIBLE_DEVICES=-1

USER appuser
ENTRYPOINT ["python", "/scripts/run.py"]
