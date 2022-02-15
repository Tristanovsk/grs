# source ~/set_proxy.sh
#
module load rclone


list_file="exe/list/list_grs_jirau.csv"
year=2021
for tile in `awk -F ',' 'NR > 1 {print $6}'  $list_file`; do
  echo rclone sync /datalake/watcal/S2-L2GRS/$tile dropbox_harmel:/satellite/S2/cnes/$tile
   rclone sync -P /datalake/watcal/S2-L2GRS/$tile/$year dropbox_harmel:/satellite/S2/cnes/$tile/$year
done

#year=2015
#tile=33VXF
#rclone sync /datalake/watcal/S2-L2GRS/$tile/$year dropbox_harmel:/satellite/S2/cnes/$tile/$year
