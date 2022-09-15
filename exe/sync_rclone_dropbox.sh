# source ~/set_proxy.sh
#

module load rclone


list_file="exe/list/list_grs_gernez_feb2022.csv"
list_file="exe/list/list_grs_cnes_obs2mod.csv"

year=2021
for tiledate in `awk -F ',' 'NR > 1 && $1 == 1 {print $6$3}'  $list_file`; do
   tile=${tiledate:0:5}
   year=${tiledate:5:4}

   year=2021
   echo rclone sync -P --include='*.nc' /datalake/watcal/S2-L2GRS/$tile/$year dropbox_harmel:/satellite/S2/cnes/$tile/$year
   #rclone sync -P --include='*.nc' /datalake/watcal/S2-L2GRS/$tile dropbox_harmel:/satellite/S2/cnes/$tile
   rclone sync -P --include='*.nc' /datalake/watcal/S2-L2GRS/$tile/$year/09 dropbox_harmel:/satellite/S2/cnes/$tile/$year/09
done

#year=2015
#tile=33VXF
#rclone sync /datalake/watcal/S2-L2GRS/$tile/$year dropbox_harmel:/satellite/S2/cnes/$tile/$year
