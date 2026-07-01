"""
Top-level CLI entry point for cellarium-cas.

Usage:
    cellarium-cas annotate [OPTIONS]
    cellarium-cas benchmark flat [OPTIONS]
    cellarium-cas benchmark ontology-aware [OPTIONS]
"""

import click

from cellarium.cas.cli.annotate import annotate_command
from cellarium.cas.cli.benchmark import benchmark_group


@click.group()
def cli() -> None:
    """Cellarium Cell Annotation Service (CAS) command-line interface."""


cli.add_command(annotate_command)
cli.add_command(benchmark_group)
