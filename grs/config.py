''' Set grs absolute path and global variables'''

import os, sys
import logging


root = os.path.dirname(os.path.abspath(__file__))

# -----------------
# please set the following path according to your specific file tree
grs_root = root

try:
    logging.info("data_root is "+os.environ.get('DATA_ROOT'))
except Exception as error:
    logging.info("environment variable DATA_ROOT does not exist")
     
data_root = os.path.abspath(os.environ.get('DATA_ROOT', '/DATA/git/satellite_app/grsdata'))
lut_root = os.path.join(data_root, "LUT")
cams_folder = os.path.join(os.environ.get('CAMS_PATH', '/datalake/watcal/ECMWF/CAMS'))  # os.path.join(data_root, "CAMS")


# -----------------
# do not change:
fortran_root = os.path.join(root, "fortran")
sys.path.insert(0, fortran_root)
