# The PTools library

A fully pythonic library for building and analysing macromolecular assemblies.

* Free software: MIT license

## Installation

We strongly recommend to install the package using [uv](https://github.com/astral-sh/uv).
This will ensure that all dependencies are correctly installed.

To install `uv`, run the following command:

```bash
# On macOS and Linux.
curl -LsSf https://astral.sh/uv/install.sh | sh
```

Then, setup ptools using the following command:

```bash
git clone https://github.com/ptools/ptools.git
cd ptools
uv sync
source .venv/bin/activate
```

PTools is ready to use!

```bash
uv run ptools --help
```


# Collaborators

- Benoist LAURENT
- Chantal PREVOST
- Hubert SANTUZ
- Charles ROBERT

Past contributors:

- Adrien SALADIN
- Benjamin BOYER
- Pierre POULAIN


# Credits

This package was created with [Cookiecutter](https://github.com/audreyr/cookiecutter)
and the [audreyr/cookiecutter-pypackage](https://github.com/audreyr/cookiecutter-pypackage)
project template.

