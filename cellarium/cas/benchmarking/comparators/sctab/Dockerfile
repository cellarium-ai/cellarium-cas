FROM python:3.10-slim

# git is required to pip-install scTab from its GitHub repo; the slim image has none.
RUN apt-get update && apt-get install -y --no-install-recommends git && rm -rf /var/lib/apt/lists/*

# Torch pinned to the 2.0.x line scTab was trained with (NVIDIA Merlin 23.02
# container ships torch 1.14.0a0 ~= torch 2.0 pre-release; scTab's own
# requirements.txt floors are torch>=1.12, lightning>2.0.0, torchmetrics>=0.11.1).
# cellnet/models.py does `import lightning.pytorch`, which is provided by the
# `lightning` distribution (pytorch-lightning alone ships no `lightning/` namespace),
# so both are installed. pytorch-lightning and lightning are pinned to 2.0.1 and
# resolved in the same command as torch, so the resolver cannot upgrade torch off the
# pinned CPU build. numpy is pinned <2 because the torch 2.0.1 wheel is built
# against NumPy 1.x and aborts at import under NumPy 2.
# --extra-index-url (not --index-url): PyPI stays the primary index, so torch's
# runtime deps (typing-extensions, jinja2, ...) resolve there; only torch==2.0.1+cpu
# comes from the CPU index (the +cpu local tag is exclusive to it).
RUN pip install --no-cache-dir --only-binary :all: \
    --extra-index-url https://download.pytorch.org/whl/cpu \
    "torch==2.0.1+cpu" \
    "pytorch-lightning==2.0.1" \
    "lightning==2.0.1" \
    "lightning-cloud==0.5.60" \
    "torchmetrics==0.11.1" \
    "numpy==1.26.4" \
    "pandas==2.3.3" \
    "scipy==1.15.3" \
    "pyarrow==25.0.1" \
    "anndata==0.11.4" \
    "pyyaml==6.0.3" \
    "click==8.5.0" \
    "tqdm==4.70.1"

# scTab is a source install pinned to commit 5ede7f2: no wheel exists, so
# --only-binary :all: cannot apply (Sonar flags this line; the pinned commit
# SHA is the version lock for it).
RUN pip install --no-cache-dir --no-deps "git+https://github.com/theislab/scTab.git@5ede7f2"

COPY run.py /scripts/run.py

RUN useradd -m -u 1000 appuser && chown -R appuser:appuser /scripts

USER appuser
ENTRYPOINT ["python", "/scripts/run.py"]