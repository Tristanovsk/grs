import os, sys
import pandas as pd
import subprocess
from multiprocessing import Pool


def call(command):
    print(command)
    pipeline_out = subprocess.call(command, stderr=subprocess.STDOUT, shell=True)
    return


ncore=8
ifile = 'tmp_grslist_220725_103950'
command = pd.read_csv(ifile).values

with Pool(processes=ncore) as pool:
    pool.map(call, command)
    pool.close
