"""
Shared filesystem path utilities.
"""

import typing as t
from pathlib import Path


def resolve_within(base: t.Union[str, Path], *parts: str) -> Path:
    """
    Join *parts* onto *base*, resolve the result, and assert it stays within *base*.

    This is the SonarCloud-recognised sanitizer pattern for path injection: after resolving
    the candidate path, :meth:`~pathlib.Path.relative_to` raises :class:`ValueError` if the
    path escaped *base* via traversal sequences (e.g. ``../../``).

    :param base: Trusted base directory.  Need not exist yet.
    :param parts: Path components to join (e.g. a plain filename).
    :returns: Resolved :class:`~pathlib.Path` guaranteed to be within *base*.
    :raises ValueError: If the resolved candidate escapes *base*.
    """
    resolved_base = Path(base).resolve()
    candidate = resolved_base.joinpath(*parts).resolve()
    try:
        candidate.relative_to(resolved_base)
    except ValueError:
        raise ValueError(f"Path '{candidate}' is outside the target directory '{resolved_base}'")
    return candidate
