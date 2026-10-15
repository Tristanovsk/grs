'''
Configuration of GRS: settings, folder of the look-up tables (LUTs) and their download from Zenodo.

The settings are read from, in increasing order of priority:

1. the default configuration shipped with the package (``grs/config.yml``),
2. the user configuration file: the file given by the ``GRS_CONFIG`` environment variable or, if not
   set, ``$XDG_CONFIG_HOME/grs/config.yml`` (``~/.config/grs/config.yml`` by default),
3. the ``GRS_DATA_ROOT`` environment variable, which overrides ``path: grsdata``.

If ``grsdata`` is not set (``null``), the LUTs are read from the default folder of the user data,
:py:func:`default_data_root` (``~/.local/share/grs/lut`` on Linux): downloading or copying the LUTs
there needs no configuration.

The user file only needs the keys that differ from the defaults, e.g.::

    path:
      grsdata: '/path/to/lut'

Create it with ``grs_config --init /path/to/lut``, download the LUTs with ``grs_config --download``
and run ``grs_config`` to print the configuration in use and check the LUTs.
'''

import argparse
import datetime
import hashlib
import json
import os
import sys
import urllib.request
from functools import lru_cache
from importlib.resources import files
from pathlib import Path

import platformdirs
import yaml

# keys of `path` giving the LUT file names, with their short name
LUT_KEYS = {'toa': 'toa_lut', 'transmittance': 'trans_lut'}


def default_data_root():
    '''
    Default folder of the LUTs, used when ``grsdata`` is not set: ``lut`` in the user data folder of
    GRS (``~/.local/share/grs/lut`` on Linux, ``~/Library/Application Support/grs/lut`` on macOS,
    ``%LOCALAPPDATA%\\grs\\lut`` on Windows).
    '''
    return platformdirs.user_data_path('grs') / 'lut'


def default_config_file():
    '''Return the packaged default configuration file.'''
    return files(__package__).joinpath('config.yml')


def user_config_file():
    '''Return the path of the user configuration file (which may not exist).'''
    if os.environ.get('GRS_CONFIG'):
        return Path(os.environ['GRS_CONFIG']).expanduser()
    config_home = os.environ.get('XDG_CONFIG_HOME') or Path.home() / '.config'
    return Path(config_home) / 'grs' / 'config.yml'


def read_user_config():
    '''
    Return the settings of the user configuration file, ``{}`` if there is none.

    :raises FileNotFoundError: if ``GRS_CONFIG`` points to a missing file
    '''
    user_file = user_config_file()
    if user_file.is_file():
        with open(user_file, 'r', encoding='utf-8') as file:
            return yaml.safe_load(file) or {}
    if os.environ.get('GRS_CONFIG'):
        raise FileNotFoundError(f'GRS_CONFIG points to a missing file: {user_file}')
    return {}


def _merge(base, update):
    '''Recursively merge the dictionary `update` into a copy of `base`.'''
    merged = dict(base)
    for key, value in update.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = _merge(merged[key], value)
        else:
            merged[key] = value
    return merged


@lru_cache(maxsize=None)
def load_config():
    '''
    Load the GRS configuration (defaults, then user file, then environment variables).

    The result is cached: call ``load_config.cache_clear()`` after changing the user file or the
    environment variables in a running session.
    '''
    with default_config_file().open('r', encoding='utf-8') as file:
        config = yaml.safe_load(file) or {}

    config = _merge(config, read_user_config())

    if os.environ.get('GRS_DATA_ROOT'):
        config = _merge(config, {'path': {'grsdata': os.environ['GRS_DATA_ROOT']}})

    grsdata = config['path'].get('grsdata')
    config['path']['grsdata'] = os.path.expanduser(grsdata) if grsdata else str(default_data_root())
    return config


def data_root_source():
    '''Return where the LUT folder in use is set: ``GRS_DATA_ROOT``, user file or default folder.'''
    if os.environ.get('GRS_DATA_ROOT'):
        return 'GRS_DATA_ROOT'
    if (read_user_config().get('path') or {}).get('grsdata'):
        return 'user configuration'
    return 'default folder'


def lut_path(name):
    '''
    Path of a LUT of the configuration in use.

    :param name: ``'toa'`` or ``'transmittance'``
    '''
    paths = load_config()['path']
    return Path(paths['grsdata'], paths[LUT_KEYS[name]])


def lut_status():
    '''Return the LUTs of the configuration in use as a list of ``(name, path, exists)``.'''
    return [(name, lut_path(name), lut_path(name).is_file()) for name in LUT_KEYS]


def check_luts(paths=None):
    '''
    Check that the LUT files exist.

    :param paths: LUT files to check (default: the LUTs of the configuration in use)
    :raises FileNotFoundError: listing the missing LUTs and how to get them
    '''
    if paths is None:
        paths = [path for _, path, _ in lut_status()]
    missing = [str(path) for path in paths if not Path(path).is_file()]
    if missing:
        raise FileNotFoundError(
            'LUT file(s) not found: ' + ', '.join(missing)
            + f'. Download them with `grs_config --download` (into {load_config()["path"]["grsdata"]}), '
            'or set the folder holding them with `grs_config --init /path/to/lut` (or path: grsdata in '
            'the user configuration file, or GRS_DATA_ROOT); check the configuration with `grs_config`.')


def init_user_config(data_root=None, force=False):
    '''
    Write the user configuration file with the LUT folder `data_root` and, commented, the default
    settings (to uncomment and change if needed).

    :param data_root: folder of the LUTs; if None, ``grsdata`` is not set and the default folder
        :py:func:`default_data_root` is used (and created)
    :param force: overwrite an existing user configuration file
    :return: path of the file written
    :raises FileExistsError: if the file exists and `force` is False
    '''
    user_file = user_config_file()
    if user_file.exists() and not force:
        raise FileExistsError(f'{user_file} already exists (use --force to overwrite it)')

    defaults = default_config_file().read_text(encoding='utf-8')
    defaults = defaults.split('\n\n', 1)[1]  # without the header of the packaged file
    commented = ''.join(('# ' + line if line.strip() else '#\n') if not line.startswith('#') else line
                        for line in defaults.splitlines(keepends=True))
    if data_root:
        data_root_lines = f"path:\n  grsdata: '{Path(data_root).expanduser().absolute()}'\n"
    else:
        data_root_lines = (f'# LUT folder: the default one, {default_data_root()}; to use another folder:\n'
                           '# path:\n'
                           "#   grsdata: '/path/to/lut'\n")
    content = (
        f'# GRS user configuration, created by `grs_config --init` on {datetime.date.today()}.\n'
        '# Only the keys that differ from the packaged defaults are needed: to change a setting, copy\n'
        '# it (with its parent keys) from the defaults below and uncomment it.\n'
        '# Print the configuration in use and check the LUTs with `grs_config`.\n'
        '\n'
        f'{data_root_lines}'
        '\n'
        '# ------------------------------------------------------------------------------------------\n'
        '# Packaged defaults (grs/config.yml)\n'
        '# ------------------------------------------------------------------------------------------\n'
        f'{commented}'
    )
    user_file.parent.mkdir(parents=True, exist_ok=True)
    user_file.write_text(content, encoding='utf-8')
    if not data_root:
        default_data_root().mkdir(parents=True, exist_ok=True)
    load_config.cache_clear()
    return user_file


def zenodo_record_url():
    '''
    URL of the Zenodo API record of the LUTs (``zenodo`` section of the configuration).

    :raises ValueError: if no record is configured
    '''
    zenodo = load_config().get('zenodo') or {}
    if not zenodo.get('record'):
        raise ValueError('no Zenodo record of the LUTs is configured (zenodo: record); download the '
                         'LUTs manually and copy them into the LUT folder (see `grs_config`)')
    return f"{zenodo.get('url', 'https://zenodo.org').rstrip('/')}/api/records/{zenodo['record']}"


def download_luts(force=False, progress=True):
    '''
    Download the LUTs of the configuration in use from Zenodo into the LUT folder.

    Each file is checked against the MD5 checksum given by Zenodo; it is written to a ``.part`` file
    renamed at the end, so that an interrupted download leaves no incomplete LUT.

    :param force: download the files that already exist
    :param progress: print the progress of the downloads
    :return: list of the paths of the LUTs
    '''
    with urllib.request.urlopen(zenodo_record_url()) as response:
        record = json.load(response)
    remote = {f['key']: f for f in record['files']}

    paths = []
    for name, path, exists in lut_status():
        paths.append(path)
        if exists and not force:
            if progress:
                print(f'LUT {name:13s} already there: {path}')
            continue
        if path.name not in remote:
            raise FileNotFoundError(f'{path.name} is not in the Zenodo record {record["links"].get("self_html", "")}')
        entry = remote[path.name]
        algorithm, expected = entry['checksum'].split(':', 1)
        path.parent.mkdir(parents=True, exist_ok=True)
        part = path.with_name(path.name + '.part')
        digest = hashlib.new(algorithm)
        done = 0
        with urllib.request.urlopen(entry['links']['self']) as response, open(part, 'wb') as out:
            for block in iter(lambda: response.read(1 << 20), b''):
                out.write(block)
                digest.update(block)
                done += len(block)
                if progress:
                    print(f'\rLUT {name:13s} downloading {path.name}: {done / 1e6:.0f}/{entry["size"] / 1e6:.0f} MB',
                          end='', flush=True)
        if progress:
            print(flush=True)
        if digest.hexdigest() != expected:
            part.unlink()
            raise IOError(f'{algorithm} checksum mismatch for {path.name}: download it again')
        part.replace(path)
    return paths


def main(argv=None):
    '''
    ``grs_config`` command: print the configuration in use and check the LUTs (exit code 1 if a LUT
    is missing), create the user configuration file with ``--init``, download the LUTs with
    ``--download``.
    '''
    parser = argparse.ArgumentParser(
        prog='grs_config',
        description='Print the GRS configuration in use and check the LUTs, create the user '
                    'configuration file, or download the LUTs from Zenodo.')
    parser.add_argument('--init', nargs='?', const='', metavar='DATA_ROOT',
                        help='create the user configuration file, with DATA_ROOT the folder of the LUTs '
                             '(default folder if not given)')
    parser.add_argument('--download', action='store_true',
                        help='download the missing LUTs from Zenodo into the LUT folder')
    parser.add_argument('--force', action='store_true',
                        help='with --init, overwrite an existing user configuration file; '
                             'with --download, download the LUTs again')
    args = parser.parse_args(argv)

    if args.init is not None:
        try:
            user_file = init_user_config(args.init or None, force=args.force)
        except FileExistsError as error:
            print(f'error: {error}', file=sys.stderr)
            return 1
        print(f'user configuration written: {user_file}')
        print('---')

    if args.download:
        try:
            download_luts(force=args.force)
        except (ValueError, OSError) as error:
            print(f'error: {error}', file=sys.stderr)
            return 1
        print('---')

    user_file = user_config_file()
    print(f'default configuration: {default_config_file()}')
    print(f'user configuration:    {user_file} ({"found" if user_file.is_file() else "not found"})')
    print(f'LUT folder:            {load_config()["path"]["grsdata"]} ({data_root_source()})')
    print('---')
    print(yaml.safe_dump(load_config(), sort_keys=False), end='')
    print('---')
    missing = False
    for name, path, exists in lut_status():
        print(f'LUT {name:13s} {"ok     " if exists else "MISSING"} {path}')
        missing |= not exists
    if missing:
        print(f'error: LUT file(s) missing; download them with `grs_config --download`, or copy them into '
              f'{load_config()["path"]["grsdata"]}, or set their folder with `grs_config --init /path/to/lut` '
              '(or path: grsdata in the user configuration file, or GRS_DATA_ROOT)', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
