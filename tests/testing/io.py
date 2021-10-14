"""testing.io - I/O testing utilities."""

from contextlib import contextmanager
import tempfile


def random_filename() -> str:
    """Returns a random file name."""
    with tempfile.NamedTemporaryFile() as tmpfile:
        return tmpfile.name


# Ignores R1732: Consider using 'with' for resource-allocating operations
# pylint: disable=R1732


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


@contextmanager
def mk_empty_file(**kwargs) -> tempfile.NamedTemporaryFile:
    """Creates a temporary empty file."""
    try:
        tmpfile = tempfile.NamedTemporaryFile(mode="wt", **kwargs)
        yield tmpfile
    finally:
        tmpfile.close()
