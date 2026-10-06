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
conda create -n grs_dev python=3.11
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
pytest grs/tests/integration_tests/
```

This is the same test suite executed by CI (job `python-tests` in
[`.github/workflows/main.yml`](.github/workflows/main.yml)), which also requires GDAL to be
installed.

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
pip install sphinx sphinx-rtd-theme myst-parser sphinx-autoapi
cd docs
sphinx-build -b html source build/html
```

## Code of conduct

Please be respectful and constructive in issues, pull requests, and discussions. For anything
sensitive, contact the maintainers directly (see [Authors](README.md#authors) in the README).

## License

By contributing, you agree that your contributions will be licensed under the
[Apache License 2.0](LICENSE), the license of this project.
