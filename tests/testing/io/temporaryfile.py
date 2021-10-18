"""Utilities to create temporary files."""

from contextlib import contextmanager
import tempfile


@contextmanager
def mk_tmp_file(content: str = "", **kwargs) -> tempfile.NamedTemporaryFile:
    """Creates a temporary file."""
    try:
        tmpfile = tempfile.NamedTemporaryFile("wt", **kwargs)
        tmpfile.write(content)
        tmpfile.flush()
        tmpfile.seek(0)
        yield tmpfile
    finally:
        tmpfile.close()


def mk_empty_file(**kwargs) -> tempfile.NamedTemporaryFile:
    """Creates a temporary empty file."""
    return mk_tmp_file(content="")


