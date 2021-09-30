"""testing.io - I/O testing utilities."""

import tempfile


def random_filename():
    """Return a random file name."""
    with tempfile.NamedTemporaryFile() as tmpfile:
        return tmpfile.name

# Ignores R1732: Consider using 'with' for resource-allocating operations
# pylint: disable=R1732
def mk_tmp_file(content="", mode="wt", **kwargs):
    """Create a temporary empty file.
    Returns:
        tempfile.NamedTemporaryFile
    """
    tmpfile = tempfile.NamedTemporaryFile(mode, **kwargs)
    tmpfile.write(content)
    tmpfile.flush()
    tmpfile.seek(0)
    return tmpfile


def mk_empty_file(**kwargs):
    return mk_tmp_file(content="", **kwargs)
