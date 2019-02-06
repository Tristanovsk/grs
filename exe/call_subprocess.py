import os, sys
import subprocess



dir = os.path.dirname(os.path.abspath(__file__))
exe = os.path.join(dir, 'grs_proc_sites.py')
site_file = os.path.join(dir, 'aeronet-oc_sites.txt')
odir = os.path.abspath('/nfs/DP/S2/L2/GRS/aeronet-oc')
site_file = os.path.join(dir, 'amazone_sites.txt')
odir = os.path.abspath('/nfs/DP/S2/L2/GRS/amazone')
site_file = os.path.join(dir, 'gernez_sites.txt')
odir = os.path.abspath('/nfs/DP/S2/L2/GRS/gernez')


while True:
    try:
        pipeline_out = subprocess.call(['python3', exe, site_file, odir], stderr=subprocess.STDOUT)
    except:
        sys.exit()
