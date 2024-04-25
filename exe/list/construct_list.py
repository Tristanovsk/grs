import os
import glob
import pandas as pd
opj = os.path.join
idir="/datalake/watcal/L8-L1-C2"
start_date='2013-01-01'
end_date='2022-09-01'

files = glob.glob(opj(idir,'*','*','*','*','*.tar'))
files = pd.DataFrame({"file":files})
files = files.file.str.split('/',expand=True).iloc[:,1:]
files.columns=['rep1','rep2','sat','tile','year','mm','dd','file']

tiles =  files.tile.unique()
list = pd.DataFrame({'process (yes if 1)':1,
              'Site Name':'_v15','start_date':start_date,'end_date':end_date,'satellite':'landsat','tile':tiles,'resolution (m)':30,'flag':1

})
list.to_csv('./exe/list/list_landsat_jegou.csv',index=False)