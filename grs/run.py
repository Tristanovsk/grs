''' Executable to process Sentinel-2 L1C images for aquatic environment

Usage:
  grs <input_file> [-o <ofile>] [--odir <odir>]  [--output param] [--resolution=res] \
   [--levname <lev>] [--no_clobber] [--allpixels]
  grs -h | --help
  grs -v | --version

Options:
  -h --help        Show this screen.
  -v --version     Show version.

  <input_file>     Input file to be processed

  -o ofile         Full (absolute or relative) path to output L2 image.
  --odir odir      Ouput directory [default: ./]
  --levname lev    Level naming used for output product [default: L2grs]
  --no_clobber     Do not process <input_file> if <output_file> already exists.
  --output param   set output unit: 'Rrs' or 'Lwn' [default: Rrs]
  --resolution=res  spatial resolution of the scene pixels
  --allpixels      force to process all pixels whatever they are masked (cloud, vegetation...) or not

'''

import os
from docopt import docopt
import logging
from .config import *
from .grs_process import process


def main():
    args = docopt(__doc__, version=__package__ + ' ' + VERSION)
    print(args)

    file = args['<input_file>']
    lev = args['--levname']
    cams_dir = args['--cams_dir']
    noclobber = args['--no_clobber']
    allpixels = args['--allpixels']
    resolution = args['--resolution']
    output = args['--output']

    ##################################
    # File naming convention
    ##################################

    outfile = args['-o']
    if outfile == None:
        basename = os.path.basename(file)
        outfile = basename.replace('L1C', lev)
        outfile = outfile.replace('.SAFE', '.nc').rstrip('/')

    odir = args['--odir']
    if odir == './':
        odir = os.getcwd()

    if not os.path.exists(odir):
        os.makedirs(odir)

    outfile = os.path.join(odir, outfile)

    if os.path.isfile(outfile) & noclobber:
        print('File ' + outfile + ' already processed; skip!')
        sys.exit()

    logging.info('call grs_process for the following paramater. File:' +
                 file + ', output file:' + outfile +
                 ', cams_dir:' + cams_dir +
                 ', resolution:' + str(resolution))
    process().execute(file, outfile, cams_dir, resolution=resolution,
                      output=output, allpixels=allpixels)
    return


if __name__ == "__main__":
    main()
