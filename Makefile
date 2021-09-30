.PHONY: clean clean-test clean-pyc clean-build docs help
.DEFAULT_GOAL := help
define BROWSER_PYSCRIPT
import os, webbrowser, sys
try:
	from urllib import pathname2url
except:
	from urllib.request import pathname2url

webbrowser.open("file://" + pathname2url(os.path.abspath(sys.argv[1])))
endef
export BROWSER_PYSCRIPT

define PRINT_HELP_PYSCRIPT
import re, sys

for line in sys.stdin:
	match = re.match(r'^([a-zA-Z_-]+):.*?## (.*)$$', line)
	if match:
		target, help = match.groups()
		print("%-20s %s" % (target, help))
endef
export PRINT_HELP_PYSCRIPT
BROWSER := python -c "$$BROWSER_PYSCRIPT"

help:
	@python -c "$$PRINT_HELP_PYSCRIPT" < $(MAKEFILE_LIST)

clean: clean-build clean-pyc clean-test ## remove all build, test, coverage and Python artifacts


clean-build: ## remove build artifacts
	rm -fr build/
	rm -fr dist/
	rm -fr .eggs/
	find . -name '*.egg-info' -exec rm -fr {} +
	find . -name '*.egg' -exec rm -rf {} +

clean-pyc: ## remove Python file artifacts
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -fr {} +

clean-test: ## remove test and coverage artifacts
	rm -fr .tox/
	rm -f .coverage
	rm -fr htmlcov/
	rm -fr .cache/

lint: lint-package lint-tests ## check style with pylint

# Ignore C0103: Variable name "c" doesn't conform to snake_case naming style
lint-package:
	pylint --disable=C0103 ptools

# Ignores:
#   - C0103: Variable name "c" doesn't conform to snake_case naming style
#   - C0114: Missing module docstring
#   - C0115: Missing class docstring
#   - C0116: Missing function or method docstring
lint-tests:
	pylint --disable=C0103,C0114,C0115,C0116 tests

test: ## run tests quickly with the default Python
	py.test --ignore=tests/pyattract/test_attract.py
	@echo "\033[33m** WARNING: didn't test tests/pyattract/test_attract.py **\033[0m"

test-all: ## run tests on every Python version with tox
	tox

coverage: ## check code coverage quickly with the default Python
	py.test --ignore=tests/pyattract/test_attract.py --cov-report term-missing --cov=ptools tests/
	@echo "\033[33m** WARNING: didn't test tests/pyattract/test_attract.py **\033[0m"


docs: ## generate Sphinx HTML documentation, including API docs
	rm -f docs/ptools.rst
	rm -f docs/modules.rst
	sphinx-apidoc -o docs/ ptools
	$(MAKE) -C docs clean
	$(MAKE) -C docs html
	$(BROWSER) docs/_build/html/index.html

servedocs: docs ## compile the docs watching for changes
	watchmedo shell-command -p '*.rst' -c '$(MAKE) -C docs html' -R -D .

release: clean ## package and upload a release
	python setup.py sdist upload
	python setup.py bdist_wheel upload

dist: clean ## builds source and wheel package
	python setup.py sdist
	python setup.py bdist_wheel
	ls -l dist

install: clean ## install the package to the active Python's site-packages
	python setup.py install
