import os
import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

plt.rcParams.update({'font.size': 16})

file='issues/timelag/S2_timelag_MSI_rel.csv'
lags = pd.read_csv(file,header=None,skiprows=1,index_col=0)
lags = lags.T.sort_values('band')
reflfile =  'issues/timelag/mega_glinted_pixels.csv'
refl = pd.read_csv(reflfile,sep='\t')
refl['lags']= lags.timelag.values.astype('float')
refl= refl.set_index(['Wavelength','lags'])

fig, axs = plt.subplots(figsize=(8, 5))
#axs = axs.ravel()
for name, refl_ in  refl.iteritems():
    print(refl_.index)
    wl =refl_.index._get_level_values(0)
    timelag= refl_.index._get_level_values(1)
    axs.plot(wl, refl_, label=name , alpha=0.75)
    c =axs.scatter(wl, refl_,  cmap="Spectral_r",c=timelag, lw=0.5, alpha=1)
plt.legend()

divider = make_axes_locatable(axs)
cax = divider.append_axes('right', size='5%', pad=0.05)
cbar = fig.colorbar(c,cax=cax, format=mpl.ticker.ScalarFormatter(),
                    shrink=1.0, fraction=0.1, pad=0)

axs.set_ylabel('TOA reflectance (L1C)')
axs.set_xlabel('Wavelength (nm)')



