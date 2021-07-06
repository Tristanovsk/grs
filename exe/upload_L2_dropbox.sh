
for tile in 30SUF 30TWT 30TXR 30TWS 30TWT 18SVH 31UET 29SPB 34HBH 17RLK 30UYV 34VDK 34VCJ 55KDV 34SFH 41RKJ 31TFJ 18GYU 17RLK; do
  dropbox_uploader.sh -s upload /datalake/watcal/S2-L2GRS/$tile /satellite/S2/cnes/
done

for tile in 35SMD 36SWD; do
  dropbox_uploader.sh -s upload /datalake/watcal/S2-L2GRS/$tile /satellite/S2/cnes/
done


