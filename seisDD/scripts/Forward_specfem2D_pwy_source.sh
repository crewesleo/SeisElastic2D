#!/bin/bash

isource=$1
NPROC_SPECFEM=$2
data_tag=$3
data_list=$4
velocity_dir=$5
SAVE_FORWARD=$6
WORKING_DIR=$7
DISK_DIR=$8
DATA_DIR=$9
job=${10}

if [ $isource -eq 1 ] ; then
    echo "SPECFEM2D Forward Modeling ..."
    echo "NPROC_SPECFEM=$NPROC_SPECFEM"
    echo "data_tag=$data_tag"
    echo "data_list=$data_list"
    echo "velocity_dir=$velocity_dir"
    echo "SAVE_FORWARD=$SAVE_FORWARD"
    echo "WORKING_DIR=$WORKING_DIR"
    echo "DISK_DIR=$DISK_DIR"
    echo "DATA_DIR=$DATA_DIR"
    echo "job=$job"
fi

ISRC_WORKING_DIR=$( seq --format="$WORKING_DIR/%06.f/" $(($isource-1)) $(($isource-1)) ) # working directory (on local nodes, where specfem runs)
ISRC_DATA_DIR=$( seq --format="$DISK_DIR/%06.f/" $(($isource-1)) $(($isource-1)) )/$data_tag

rm -rf $ISRC_WORKING_DIR $ISRC_DATA_DIR
mkdir -p $ISRC_WORKING_DIR $ISRC_DATA_DIR

cd $ISRC_WORKING_DIR

##echo "####### copy executables & input files ######"
cp -r $SUBMIT_DIR/parameter ./
cp -r $SUBMIT_DIR/bin ./
cp -r $SUBMIT_DIR/DATA ./

## copy STATIONS
ISRC_WORKING_DIR_STATION=$ISRC_WORKING_DIR/DATA/SOURCES
ISRC_WORKING_DIR_STATION_isource=$( seq --format="$ISRC_WORKING_DIR_STATION/%06.f/" $(($isource-1)) $(($isource-1)) ) # STATIONS corresponding to current isource

cp -r $ISRC_WORKING_DIR_STATION_isource/*.dat ./DATA

# if $SUBMIT_DIR/SU_process exist
if [ -d "$SUBMIT_DIR/SU_process" ]; then
#    cp -r $SUBMIT_DIR/SU_process ./
    cp -r $ISRC_WORKING_DIR_STATION_isource/SU_process ./
fi
# if velocity_dir exist and not empty (PWY)
#if [ "$(ls -A $velocity_dir)" ]; then
echo "try velocity directory: $velocity_dir"
cp -r $velocity_dir/* ./DATA
#fi

# PWY
#rm -rf OUTPUT_FILES SEM
mkdir -p  OUTPUT_FILES SEM
rm -rf ./OUTPUT_FILES/*.su
# Source location
export xs=$(awk -v "line=$isource" 'NR==line { print $1 }' DATA/sources.dat)
export zs=$(awk -v "line=$isource" 'NR==line { print $2 }' DATA/sources.dat)

### echo " edit SOURCE "
FILE="./DATA/SOURCE"
sed -e "s/^xs.*$/xs =    $xs/g" $FILE > temp;  mv temp $FILE
sed -e "s/^zs.*$/zs =    $zs/g" $FILE > temp;  mv temp $FILE

##### edit 'Par_file' #####
FILE="./DATA/Par_file"
sed -e "s#^SIMULATION_TYPE.*#SIMULATION_TYPE = 1 #g"  $FILE > temp; mv temp $FILE
sed -e "s#^SAVE_FORWARD.*#SAVE_FORWARD = .$SAVE_FORWARD. #g"  $FILE > temp; mv temp $FILE


##### forward simulation (data) #####
./bin/xmeshfem2D

if [ $isource -eq 1 ] ; then
    echo "mpirun -np $NPROC_SPECFEM ./bin/xspecfem2D"
fi
mpirun -np $NPROC_SPECFEM ./bin/xspecfem2D

## copy and preprocessing of data 
arr=$(echo $data_list | tr "," "\n")

for x in $arr
do
    if [ -f "SU_process/process_syn.sh" ]; then
        sh SU_process/process_syn.sh \
            OUTPUT_FILES/U${x}_file_single.su \
            $ISRC_DATA_DIR/U${x}_file_single.su
    else
        cp OUTPUT_FILES/U${x}_file_single.su \
            $ISRC_DATA_DIR/U${x}_file_single.su
    fi    
done

if [ "$data_tag" == "DATA_obs" ] && [ "$job" == "modeling" ]; 
then
    mkdir -p $DATA_DIR 
    ISRC_DATA_DIR_SAVE=$( seq --format="$DATA_DIR/%06.f/" $(($isource-1)) $(($isource-1)) )
    rm -rf $ISRC_DATA_DIR_SAVE
    mkdir -p $ISRC_DATA_DIR_SAVE 
    cp -r $ISRC_DATA_DIR/* $ISRC_DATA_DIR_SAVE/
fi
