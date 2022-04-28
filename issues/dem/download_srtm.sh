

odir="/home/eh/harmelt/.snap/auxdata/dem/SRTM\ 3Sec"
eval cd $odir

for  i in $(seq -f "%02g" 1 99);do #((i=1;i<100;i++));do
  for j in $(seq -f "%02g" 1 99);do # ((j=1;j<100;j++));do
   echo wget https://download.esa.int/step/auxdata/dem/SRTM90/tiff/srtm_"$i"_"$j".zip
   wget https://download.esa.int/step/auxdata/dem/SRTM90/tiff/srtm_"$i"_"$j".zip
  done
done