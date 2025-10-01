# The PTools library

A fully pythonic library for building and analysing macromolecular assemblies.

* Free software: MIT license

## Installation

### Installing `uv`

We strongly recommend to install the package using [uv](https://github.com/astral-sh/uv).
This will ensure that all dependencies are correctly installed.

To install `uv`, run the following command:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

More details available at [uv](https://github.com/astral-sh/uv).


### Installing PTools

To install PTools, run the following command:

```bash
uv venv ptools  # setup a virtual environment named 'ptools'
source ptools/.venv/bin/activate  # activate the virtual environment
uv pip install git+https://github.com/ptools/ptools.git  # install PTools
```

### Setting up PTools development environment

If you want to contribute to PTools, you can setup a development environment using the following command:

```bash
git clone https://github.com/ptools/ptools.git
cd ptools
uv sync
```

PTools is ready to use!

```bash
uv run ptools --help
```

If you're using the library, you can run your script using the following command:

```bash
uv run my_script.py
```

To work outside of PTools source directory, you need to activate the virtual environment:

```bash
source /path/to/ptools/.venv/bin/activate
```


## Contributors

- Benoist LAURENT
- Chantal PREVOST
- Hubert SANTUZ
- Charles ROBERT

Past contributors:

- Adrien SALADIN
- Benjamin BOYER
- Pierre POULAIN
