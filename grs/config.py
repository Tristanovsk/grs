''' Set grs absolute path and global variables'''

import os, sys
from .__init__ import __version__ as VERSION
from .__init__ import __package__

root = os.path.dirname(os.path.abspath(__file__))

# -----------------
# please set the following path according to your specific file tree
# grs_root = os.path.abspath('/DATA/S2_processing/PYTHON/grs/grs')
grs_root = root

# data_root = os.path.abspath('/DATA/S2_processing/PYTHON/grs/grs/../..')
data_root = os.path.abspath('/home/grs2/grsdata')
# directory to store temporary unzipped files
tmp_dir = os.path.abspath('/tmp')
lut_root = os.path.join(data_root, "LUT")
cams_folder = os.path.abspath('/home/CAMS')  # os.path.join(data_root, "CAMS")
smac_root = os.path.join(data_root, "SMAC_COEFS")

# -----------------
# do not change:
fortran_root = os.path.join(root, "fortran")
sys.path.insert(0, fortran_root)
