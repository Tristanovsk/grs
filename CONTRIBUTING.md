# Contributing to GRSProcessor

Thanks for your interest in contributing! This document explains how to report issues, set up
a development environment, and submit changes.

## Reporting issues

Please use [GitHub Issues](https://github.com/CNES/GRSprocessor/issues) to report bugs or
request features. Include:

- the `grs` version (`grs -v`) and how you installed it (pip, conda, Docker),
- the exact command line used,
- the relevant excerpt of `log_file.log` / `error.log` (see the
  [processing chain documentation](https://cnes.github.io/GRSprocessor/processing_chain.html#log-format)
  for the log format),
- if possible, a minimal way to reproduce (input product type/sensor, resolution...).

## Development setup

```bash
git clone https://github.com/CNES/GRSprocessor.git
cd GRSprocessor
conda create -n grs_dev python=3.12
conda activate grs_dev
pip install -r requirements.txt
pip install -e .[dev]
```

The `dev` extra (see `pyproject.toml`) installs `black`, `isort`, `bumpver`, `pip-tools` and
`pytest`.

See the [README](README.md) for details on the LUT data (`grsdata`) and the `config.yml` file
required to actually run `grs` on data.

## Code style

- Format code with `black` and `isort` before committing.
- CI also runs `ruff` and `mypy` on pull requests (see
  [`.github/workflows/main.yml`](.github/workflows/main.yml), job `lint`); please fix warnings
  they raise on the lines you touch.

## Running tests

```bash
pytest grs/tests/unit/
```

The unit tests (`grs/tests/unit/`, about 10 s) run on a small L1C subset stored in the repository
(`grs/tests/data`, made with `grs/tests/data/make_test_data.py`); their end-to-end part needs the
look-up tables and is skipped without them: set `GRS_TEST_GRSDATA` to the folder of the LUTs if it is not
the `grsdata` folder of `grs/config.yml`.

This is the same test suite executed by CI (job `python-tests` in
[`.github/workflows/main.yml`](.github/workflows/main.yml)).

## Submitting changes

This project follows [Gitflow](https://nvie.com/posts/a-successful-git-branching-model/) branch
naming:

- `feature/<short-description>` — off `develop`, merged back into `develop`. Use this for most
  changes.
- `bugfix/<short-description>` — off `develop`, merged back into `develop`. Use this for
  non-urgent bug fixes.
- `release/<version>` — off `develop`, merged into **both** `main` and `develop` once ready.
- `hotfix/<version>` — off `main`, merged into **both** `main` and `develop` once ready. Use this
  for urgent fixes to a released version.

Steps:

1. Create a branch off `develop` (or off `main` for a `hotfix/*`) using one of the prefixes above.
2. Make your changes, with tests where relevant.
3. Open a pull request targeting `develop` (`main` too for `release/*`/`hotfix/*`, per Gitflow —
   open a second PR). CI runs automatically on PRs to `main` and `develop`. Make sure the
   `python-tests` job passes.
4. One of the maintainers will review your PR.

## Documentation

The documentation lives under [`docs/`](docs/) (Sphinx) and is published to
[cnes.github.io/GRSprocessor](https://cnes.github.io/GRSprocessor/) on every tag push.
See `docs/source/index.rst` for the entry point. You can build it locally with:

```bash
pip install ".[docs]"
sphinx-build -b html docs/source docs/build/html
```

## Releasing on PyPI

GRS is published on PyPI as [`GRSprocessor`](https://pypi.org/project/GRSprocessor/) (the name `grs` is
taken by another project); the package and the commands keep the name `grs`. The LUTs are not in the
package: they are downloaded from [Zenodo](https://doi.org/10.5281/zenodo.23394202) on the first run
(`grs.config.ensure_luts`, or `grs_config --download`). A new LUT version means a new Zenodo record: update
`zenodo: record` and the LUT file names in `grs/config.yml`.

1. Set the version in `grs/__init__.py` (`__version__`, the only place: `pyproject.toml` reads it) and
   move the `[Unreleased]` entries of `CHANGELOG.md` under the new version, with the date.
2. Build from a clean checkout (old `build/` or `*.egg-info` folders would bring stale files into the
   package), and check the result:

   ```bash
   git clone --branch <branch> https://github.com/CNES/GRSprocessor.git /tmp/grs_release
   cd /tmp/grs_release
   python -m pip install build twine
   python -m build                 # dist/grsprocessor-<version>.tar.gz and -py3-none-any.whl
   python -m twine check dist/*
   unzip -l dist/*.whl             # about 2 MB: code, config.yml, gases LUTs, auxiliary data; no tests
   ```

3. Try the wheel in a new environment, LUT download included:

   ```bash
   python -m venv /tmp/grs_venv && . /tmp/grs_venv/bin/activate
   pip install dist/grsprocessor-<version>-py3-none-any.whl
   GRS_DATA_ROOT=/tmp/grs_lut grs_config --download && GRS_DATA_ROOT=/tmp/grs_lut grs_config
   grs -v
   ```

4. Upload to TestPyPI first, check the project page and an installation from it, then to PyPI (API
   token of the `GRSprocessor` project in `~/.pypirc`, or `TWINE_USERNAME=__token__` and
   `TWINE_PASSWORD=<token>`):

   ```bash
   python -m twine upload --repository testpypi dist/*
   pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ GRSprocessor==<version>
   python -m twine upload dist/*
   ```

5. Tag the release (`git tag v<version> && git push --tags`).

The workflow [`.github/workflows/pypi.yml`](.github/workflows/pypi.yml) does steps 2 and 4 on GitHub
when a `v*` tag is pushed (or by hand from the Actions tab), with PyPI trusted publishing: no token is
stored, but the GitHub repository and the workflow must first be declared as a trusted publisher of the
`GRSprocessor` project on pypi.org (and on test.pypi.org).

## Code of conduct

Please be respectful and constructive in issues, pull requests, and discussions. For anything
sensitive, contact the maintainers directly (see [Authors](README.md#authors) in the README).

## License

By contributing, you agree that your contributions will be licensed under the
[Apache License 2.0](LICENSE), the license of this project.
