"""Command line (grs.run): options and defaults."""

from docopt import docopt

from grs import run


def parse(*args):
    return docopt(run.__doc__, argv=list(args))


def test_defaults():
    args = parse('S2A_MSIL1C_20210211T104151_N0500_R008_T31TFJ_20230608T124831.SAFE', '--cams_file', 'cams.nc')
    assert args['--resolution'] == '60'
    assert args['--levname'] == 'L2AGRS'
    assert args['--max_cloud_cover'] == '1'
    assert args['--scale_aot'] == '1'
    assert args['--odir'] is None
    assert not args['--snap_compliant']


def test_odir_alias():
    for option in ('-o', '--odir'):
        assert parse('image.SAFE', option, '/out')['--odir'] == '/out'


def test_options():
    args = parse('image.SAFE', '--resolution', '20', '--opac_model', 'DESE_rh70', '--no_clobber', '--allpixels',
                 '--dem_file', 'dem.tif', '--suffix', '_test')
    assert args['--resolution'] == '20'
    assert args['--opac_model'] == 'DESE_rh70'
    assert args['--no_clobber'] and args['--allpixels']
    assert (args['--dem_file'], args['--suffix']) == ('dem.tif', '_test')
