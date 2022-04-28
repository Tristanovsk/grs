# source ~/set_proxy.sh
#

module load rclone


list_file="exe/list/list_grs_gernez_feb2022.csv"
#list_file="exe/list/list_grs_cnes_seine.csv"

year=2022
for tiledate in `awk -F ',' 'NR > 1 && $1 == 1 {print $6$3}'  $list_file`; do
   tile=${tiledate:0:5}
   year=${tiledate:5:4}
   echo rclone sync -P --include='*.nc' /datalake/watcal/S2-L2GRS/$tile/$year dropbox_harmel:/satellite/S2/cnes/$tile/$year
   #rclone sync -P --include='*.nc' /datalake/watcal/S2-L2GRS/$tile dropbox_harmel:/satellite/S2/cnes/$tile
   rclone sync -P --include='*.nc' /datalake/watcal/S2-L2GRS/$tile/$year dropbox_harmel:/satellite/S2/cnes/$tile/$year
done

#year=2015
#tile=33VXF
#rclone sync /datalake/watcal/S2-L2GRS/$tile/$year dropbox_harmel:/satellite/S2/cnes/$tile/$year
