#!/bin/bash
echo The ID of the Docker image is........"$1"
echo The image to be processed is........."$2"
echo The CAMS files are located at........"$3"
echo The output file will be named........"$4"
echo The output file will located at......"$5"
echo The requested resolution is.........."$6"
echo The SURFWATER files are located at..."$7"

# Create and parametrize container
sudo docker run -d --name grs2 $1
sudo docker cp $2 grs2:/home/L1C
sudo docker cp $3 grs2:/home/CAMS

# Launch processing and retreive result
if [ -z "$7" ]
then
    sudo docker exec grs2 grs /home/L1C --cams_file /home/CAMS -o $4 --resolution $6
else
    sudo docker exec grs2 grs /home/L1C --surfwater $7 --cams_file /home/CAMS -o $4 --resolution $6
fi
sudo docker cp grs2:/home/$4 $5/$4
echo Result written at location $5/$4

# Delete container
sudo docker kill grs2
sudo docker rm grs2
