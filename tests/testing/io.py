"""testing.io - I/O testing utilities."""

import tempfile


def random_filename():
    """Return a random file name."""
    tmpfile = tempfile.NamedTemporaryFile()
    tmpfile.close()
    return tmpfile.name


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
