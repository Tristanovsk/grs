''' Set absolute path '''

import os, sys
from .__init__ import __version__ as VERSION
from .__init__ import __package__

root = os.path.dirname(os.path.abspath(__file__))

# set gpt from snap (usually in your snap folder e.g., path_to_snap/bin/gpt)
gpt = '/local/AIX/tristan.harmel/snap/bin/gpt'

# -----------------
# please set the following path according to your specific file tree
# grs_root = os.path.abspath('/DATA/S2_processing/PYTHON/grs/grs')
grs_root = root

# data_root = os.path.abspath('/DATA/S2_processing/PYTHON/grs/grs/../..')
data_root = os.path.abspath('/nfs/DD/grsdata')
# directory to store temporary unzipped files
tmp_dir = os.path.abspath('/nfs/DD/grsdata/tmp')
lut_root = os.path.join(data_root, "LUT")
cams_folder = os.path.abspath('/nfs/DP/ECMWF/CAMS')  # os.path.join(data_root, "CAMS")
smac_root = os.path.join(data_root, "SMAC_COEFS")

# -----------------
# do not change:
fortran_root = os.path.join(root, "fortran")
sys.path.insert(0, fortran_root)