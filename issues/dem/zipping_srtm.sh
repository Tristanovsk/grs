
dir=/datalake/static_aux/MNT/SRTM_90m
odir="/home/eh/harmelt/.snap/auxdata/dem/SRTM\ 3Sec"
for f in $dir/*hdr; do
  b=${f/.hdr/}
  b_=${b##*/}
  echo zip $odir/$b_.zip $b.hdr $b.tfw $b.tif
  eval zip $odir/$b_.zip $b.hdr $b.tfw $b.tif
done
