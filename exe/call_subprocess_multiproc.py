import os, sys
import pandas as pd
import subprocess
from multiprocessing import Pool


def call(command):
    print(command)
    pipeline_out = subprocess.call(command, stderr=subprocess.STDOUT, shell=True)
    return


ncore=7
ifile = '/local/AIX/nathalie.reynaud/Documents/teledec/'+\
        'processing/grs/batch_grs_ldst/batch_grs_ldst.sh'
command = pd.read_csv(ifile).values

with Pool(processes=ncore) as pool:
    pool.map(call, command)
    pool.close
