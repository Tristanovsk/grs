
for tile in 30SUF 30TWT 30TXR 30TWS 30TWT 18SVH 31UET 29SPB 34HBH 17RLK 30UYV 34VDK 34VCJ 55KDV 34SFH 41RKJ 31TFJ 18GYU 17RLK; do
  dropbox_uploader.sh -s upload /datalake/watcal/S2-L2GRS/$tile /satellite/S2/cnes/
done


list_file="exe/list_grs_cnes_gernez.csv"

list_file="exe/list_grs_gernez_juillet_2021.csv"
list_file="exe/list_grs_cnes_hafeez.csv"

for tile in `awk -F ',' '{print $6}'  $list_file`; do
  dropbox_uploader.sh -s upload /datalake/watcal/S2-L2GRS/$tile /satellite/S2/cnes/
done

year=2021
tile=31TFJ
dropbox_uploader.sh  -s -x .incomplete upload /datalake/watcal/S2-L2GRS/$tile/$year /satellite/S2/cnes/$tile
