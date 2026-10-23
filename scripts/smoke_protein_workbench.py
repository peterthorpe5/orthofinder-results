"""Run a small real DIAMOND, MAFFT and FastTree workflow without remote submissions."""

from __future__ import annotations

import argparse
import logging
import shutil
from dataclasses import replace
from pathlib import Path
from typing import Sequence

from orthofinder_interrogation_app.research.common import ResearchError, load_config
from orthofinder_interrogation_app.research.jobs import artifact_path, get_job, submit_job
from orthofinder_interrogation_app.research.runner import run_next_job

LOGGER = logging.getLogger(__name__)
UBIQUITIN = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
THIOREDOXIN = (
    "MSDKIIHLTDDSFDTDVLKADGAILVDFWAEWCGPCKMIAPILDEIADEYQGKLTVAKLNIDQNPGTAPKYGIR"
    "GIPTLLLFKNGEVAATKVGALSKGQLKEFLDANLA"
)


def main(*, argv: Sequence[str] | None = None) -> int:
    """Execute four synthetic proteins through real audited clustering and phylogeny.

    Args:
        argv: Named --config and --output command arguments.

    Returns:
        Zero for completed analyses, or two for missing tools or failed validation.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO)
    try:
        output = arguments.output.expanduser().resolve()
        if output.exists() and (not output.is_dir() or any(output.iterdir())):
            raise ResearchError("Smoke output must be a new or empty directory")
        config = replace(load_config(path=arguments.config), workspace=output)
        for name in ("diamond", "mafft", "fasttree"):
            if not shutil.which(cmd=getattr(config, name)):
                raise ResearchError(f"Missing executable: {name}")
        fasta = (
            f">ub_a\n{UBIQUITIN}\n>ub_b\n{UBIQUITIN[:-1]}A\n"
            f">trx_a\n{THIOREDOXIN}\n>trx_b\n{THIOREDOXIN[:-1]}G\n"
        )
        clustered = submit_job(
            workspace=output,
            operation="subcluster",
            fasta=fasta,
            parameters={"identity": 70, "mutual_coverage": 80, "evalue": 1e-10},
        )
        run_next_job(config=config, job_id=clustered)
        if get_job(workspace=output, job_id=clustered)["status"] != "COMPLETE":
            raise ResearchError("Real subclustering failed; inspect the worker log")
        aligned = submit_job(
            workspace=output,
            operation="align",
            parameters={
                "source_job": clustered,
                "source_artifact": "sequences.fasta",
                "reference_id": "ub_a",
                "infer_tree": True,
            },
        )
        run_next_job(config=config, job_id=aligned)
        if get_job(workspace=output, job_id=aligned)["status"] != "COMPLETE":
            raise ResearchError("Real alignment/tree inference failed; inspect the worker log")
        for artifact in ("alignment.fasta", "tree.nwk", "auspice.json", "alignment_summary.json"):
            artifact_path(workspace=output, job_id=aligned, artifact=artifact)
        LOGGER.info(
            "Real scientific smoke completed: subclusters=%s alignment=%s", clustered, aligned
        )
        return 0
    except (ResearchError, OSError, ValueError) as error:
        LOGGER.error("Scientific smoke failed: %s", error)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
