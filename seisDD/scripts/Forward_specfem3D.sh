#!/bin/bash

isource=$1
NPROC_SPECFEM=$2
data_type=$3
velocity_dir=$4
SAVE_FORWARD=$5
WORKING_DIR=$6
DISK_DIR=$7
DATA_DIR=$8

if [ $isource -eq 1 ] ; then
    echo "SPECFEM3D Forward Modeling ..."
    echo "NPROC_SPECFEM=$NPROC_SPECFEM"
    echo "data_type=$data_type"
    echo "velocity_dir=$velocity_dir"
    echo "SAVE_FORWARD=$SAVE_FORWARD"
    echo "WORKING_DIR=$WORKING_DIR"
    echo "DISK_DIR=$DISK_DIR"
    echo "DATA_DIR=$DATA_DIR"
fi

ISRC_WORKING_DIR=$( seq --format="$WORKING_DIR/%06.f/" $(($isource-1)) $(($isource-1)) ) # working directory (on local nodes, where specfem runs)
ISRC_DATA_DIR=$( seq --format="$DISK_DIR/%06.f/" $(($isource-1)) $(($isource-1)) )/$data_type

mkdir -p $ISRC_WORKING_DIR $ISRC_DATA_DIR

cd $ISRC_WORKING_DIR

mkdir -p  OUTPUT_FILES OUTPUT_FILES/DATABASES_MPI SEM

##echo "####### copy executables & input files ######"
cp -r $SUBMIT_DIR/bin ./
cp -r $SUBMIT_DIR/DATA ./
if [ "$(ls -A $velocity_dir)" ]; then
    cp $velocity_dir/* OUTPUT_FILES/DATABASES_MPI/
fi

# Source location
export lat=$(awk -v "line=$isource" 'NR==line { print $1 }' DATA/sources.dat)
export lon=$(awk -v "line=$isource" 'NR==line { print $2 }' DATA/sources.dat)
export dep=$(awk -v "line=$isource" 'NR==line { print $3 }' DATA/sources.dat)

if $DISPLAY_DETAILS ; then
    echo "source $isource -- location lat=$lat m lon=$lon m depth=$dep km "
fi

##### edit 'Par_file' #####
FILE="./DATA/Par_file"
sed -e "s#^SIMULATION_TYPE.*#SIMULATION_TYPE = 1 #g"  $FILE > temp; mv temp $FILE
sed -e "s#^SAVE_FORWARD.*#SAVE_FORWARD = .$SAVE_FORWARD. #g"  $FILE > temp; mv temp $FILE

if [ `grep ^USE_FORCE_POINT_SOURCE   ./DATA/Par_file | cut -d = -f 2 ` ]; then
    FILE="DATA/FORCESOLUTION"
else
    FILE="DATA/CMTSOLUTION"
fi
echo "edit SOURCE : $FILE "
sed -e "s/^latorUTM.*$/latorUTM:    $lat/g" $FILE > temp;  mv temp $FILE
sed -e "s/^longorUTM.*$/longorUTM:    $lon/g" $FILE > temp;  mv temp $FILE
sed -e "s/^depth.*$/depth:    $dep/g" $FILE > temp;  mv temp $FILE

##### forward simulation (data) #####
# creates and decomposes mesh
if [ $isource -eq 1 ] ; then
    echo "running mesher..."
fi
mpirun -np $NPROC_SPECFEM ./bin/xmeshfem3D

# runs database generation
if [ $isource -eq 1 ] ; then
    echo "running database generation..."
fi
mpirun -np $NPROC_SPECFEM ./bin/xgenerate_databases

# runs simulation
if [ $isource -eq 1 ] ; then
    echo "running solver..."
fi
mpirun -np $NPROC_SPECFEM ./bin/xspecfem3D
mv  OUTPUT_FILES/output_solver.txt  OUTPUT_FILES/output_forward.txt

# save
cp OUTPUT_FILES/*SU    $ISRC_DATA_DIR/
if [ "$data_type" == "DATA_obs" ];
then
    mkdir -p $DATA_DIR
    ISRC_DATA_DIR_SAVE=$( seq --format="$DATA_DIR/%06.f/" $(($isource-1)) $(($isource-1)) )
    rm -rf $ISRC_DATA_DIR_SAVE
    mkdir -p $ISRC_DATA_DIR_SAVE
    cp -r $ISRC_DATA_DIR/* $ISRC_DATA_DIR_SAVE/
fi
