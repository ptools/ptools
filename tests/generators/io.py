from contextlib import contextmanager
import tempfile


def generate_random_filename() -> str:
    """Returns a random file name."""
    with tempfile.NamedTemporaryFile() as tmpfile:
        return tmpfile.name


@contextmanager
def generate_tmp_file(content: str = "", **kwargs) -> tempfile.NamedTemporaryFile:
    """Creates a temporary file."""
    try:
        tmpfile = tempfile.NamedTemporaryFile("wt", **kwargs)
        tmpfile.write(content)
        tmpfile.flush()
        tmpfile.seek(0)
        yield tmpfile
    finally:
        tmpfile.close()


def generate_empty_file(**kwargs) -> tempfile.NamedTemporaryFile:
    """Creates a temporary empty file."""
    return generate_tmp_file(content="", **kwargs)
