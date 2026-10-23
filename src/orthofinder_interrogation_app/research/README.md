# Protein research engine

Generic workflow engine adapted from the supplied E3 app v0.19.1. It depends on
this package's viewer assets and export helpers, without E3 resource assumptions.
The OrthoFinder-specific membership adapter lives in `workbench_adapter.py` in
the parent package. Scientific jobs use an external writable workspace.

Read `docs/protein_workbench.md` and `docs/deepclust_exchange.md` from the repository
root. Entry point: `orthofinder-protein-workbench`. Tests: `tests/research/` and
the root workbench integration/smoke tests.
