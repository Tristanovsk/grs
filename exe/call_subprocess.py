import os, sys
import subprocess

dir = os.path.dirname(os.path.abspath(__file__))
exe = os.path.join(dir, 'proc_aeronet_sites.py')

while True:
    try:
        pipeline_out = subprocess.call(['python3 '+ exe], shell=True, stderr=subprocess.STDOUT)
    except:
        sys.exit()
