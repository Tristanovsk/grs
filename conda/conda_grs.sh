#!/bin/sh

# Après installation de snap, faire un lien symbolique de esasnappy vers snappy (qui se trouve dans le .snap du home

#desactivate conda s'il est activé
conda deactivate 2>/dev/null

#On récupère la connexion au proxy
#[ "$_" = "$0" ] && echo "Ce script doit être sourcé!" && exit 127

#read -s -p "Please enter your proxy password ? " _passwd || exit 127
#echo

#_user=${_user:-$USER}

#if [ -z ${_passwd} ]
#then
#	echo "No password"
#else
#
#	export http_proxy="http://${_user}:${_passwd}@proxy-surf.loc.cnes.fr:8050"
#	export https_proxy="http://${_user}:${_passwd}@proxy-surf.loc.cnes.fr:8050"

#echo 'http_proxy variable has to be set : '$http_proxy

#	export ftp_proxy="${http_proxy}"
#	export no_proxy=cnes.fr,sis.cnes.fr,gitlab.cnes.fr

#	unset _passwd
#	unset _user
#fi

#load avail module for conda, snap and jdk
module load snap/8.0 conda jdk/1.8.0_112

#switch to gcc/4.8.5
#module unload gcc

#set the variable environment for je jpy and snap use
export JDK_HOME=$JDKHOME
export JAVAHOME=$JAVA_HOME
export SNAP_HOME=$SNAPHOME
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:.


#activate conda environement python 3.8
conda activate /work/ALT/swot/aval/OBS2CO/conda_env/grs_py3.8
echo $USER

#configuration for snap
# go to the $SNAPHOME environment and configure snappy
pythonpath=$(which python)
echo $pythonpath

pippath=$(which pip)
echo $pippath

#Install snappy
#$SNAPHOME/bin/snappy-conf $pythonpath

if [ ! -d $CONDA_PREFIX/lib/python3.6/site-packages/esasnappy ]
then
	ln -s /work/ALT/swot/aval/OBS2CO/snap/snappy_3.8 $CONDA_PREFIX/lib/python3.8/site-packages/esasnappy
fi


echo $'snap.versionCheck.interval=NEVER\nsnap.jai.tileCacheSize=1024' > $HOME/.snap/etc/snap.properties

#Install the jpy package used by snappy
#pip install /work/ALT/swot/aval/OBS2CO/snap/snappy/lib/jpy-0.9.0-cp36-cp36m-linux_x86_64.whl nco

#Modifier la config jpy pour java
#cat $CONDA_PREFIX/lib/python3.6/site-packages/jpyconfig.py | grep '#java_home'

#if [ $? == 1 ]
#then
#	sed -i 's/java_home/#java_home/g' $CONDA_PREFIX/lib/python3.6/site-packages/jpyconfig.py
#	sed -i 's/jvm_dll/#jvm_dll/g' $CONDA_PREFIX/lib/python3.6/site-packages/jpyconfig.py
#fi

#if it is not already installed
#conda install GDAL
#compilation

cd /work/ALT/swot/aval/OBS2CO/git/grs2

#catch egm96 data
#mkdir -p $HOME/.snap/auxdata/dem/egm96
#cp /work/ALT/swot/aval/OBS2CO/git/grsdata/dem/ww15mgh_b.zip $HOME/.snap/auxdata/dem/egm96/ww15mgh_b.zip


if [ $# != 0 ]
then
	if [ $1 = "-c" ]
	then
        	echo "Compilation"
        	make
	else 
		if [ $1 = "-i" ]
		then
        		echo "install"
        		python setup.py install
		else
			if [ $1 = "-ci" ]
			then
				echo 'compilation and installation'
				make
				python setup.py install
			fi
		fi
	fi
fi


echo "Tester GRS sur une tuile sentinel 2 : grs /work/ALT/swot/aval/OBS2CO/git/grsdata/INPUT_DATA/S2B_MSIL1C_20180927T103019_N0206_R108_T31TGK_20180927T143835.SAFE --shape test/data/shape/SPO04.shp --odir test/results/ --aerosol cams_forecast --dem --resolution 20 \nPreciser le répertoire de sortie dans --odir"
