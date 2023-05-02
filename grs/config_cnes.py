''' Set grs absolute path and global variables'''

import os, sys
import logging
from .__init__ import __version__ as VERSION
from .__init__ import __package__

root = os.path.dirname(os.path.abspath(__file__))

# -----------------
# please set the following path according to your specific file tree
grs_root = root

try:
    logging.info("data_root is "+os.environ.get('DATA_ROOT'))
except Exception as error:
    logging.info("environment variable DATA_ROOT does not exist")
     
data_root = os.path.abspath(os.environ.get('DATA_ROOT', '/datalake/watcal/GRS/grsdata'))
lut_root = os.path.join(data_root, "LUT")

# -----------------
# do not change:
fortran_root = os.path.join(root, "fortran")
sys.path.insert(0, fortran_root)
