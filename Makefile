PYTHON ?= python3
PROJECT_ROOT := $(CURDIR)
export PYTHONPATH := $(PROJECT_ROOT)/src
export MPLCONFIGDIR := $(PROJECT_ROOT)/.matplotlib

.PHONY: test predata notebook validate all analysis

test:
	$(PYTHON) -m unittest discover -s tests -v

predata:
	$(PYTHON) -m adipose_senescence_ifn.cli build-signatures \
		--cache data/external/signatures \
		--private-out results/private \
		--public-out results/predata

notebook: predata
	jupyter nbconvert --to notebook --execute --inplace \
		--ExecutePreprocessor.timeout=300 notebooks/00_signature_overlap.ipynb

validate:
	$(PYTHON) -m adipose_senescence_ifn.cli validate-repository --root .

all: test predata notebook validate

# This target intentionally exits non-zero until matched study data pass Gate 0.
analysis:
	$(PYTHON) -m adipose_senescence_ifn.cli analyze \
		--metadata data/raw/metadata.tsv \
		--expression data/raw/expression.tsv.gz \
		--signatures data/external/signatures/mouse_signatures.json \
		--out results/private/analysis
