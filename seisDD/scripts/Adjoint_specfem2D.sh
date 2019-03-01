#!/bin/bash

isource=$1
NPROC_SPECFEM=$2
data_type=$3
velocity_dir=$4
SAVE_FORWARD=$5
WORKING_DIR=$6
DISK_DIR=$7

if [ $isource -eq 1 ] ; then
    echo "SPECFEM2D Adjoint Modeling ..."
    echo "NPROC_SPECFEM=$NPROC_SPECFEM"
    echo "data_type=$data_type"
    echo "velocity_dir=$velocity_dir"
    echo "SAVE_FORWARD=$SAVE_FORWARD"
    echo "WORKING_DIR=$WORKING_DIR"
    echo "DISK_DIR=$DISK_DIR"
fi


ISRC_WORKING_DIR=$( seq --format="$WORKING_DIR/%06.f/" $(($isource-1)) $(($isource-1)) ) # working directory (on local nodes, where specfem runs)
ISRC_DATA_DIR=$( seq --format="$DISK_DIR/%06.f/" $(($isource-1)) $(($isource-1)) )/$data_type

mkdir -p $ISRC_WORKING_DIR $ISRC_DATA_DIR

cd $ISRC_WORKING_DIR

##### edit 'Par_file' #####
FILE="./DATA/Par_file"
sed -e "s#^SIMULATION_TYPE.*#SIMULATION_TYPE = 3 #g"  $FILE > temp; mv temp $FILE
sed -e "s#^SAVE_FORWARD.*#SAVE_FORWARD = .$SAVE_FORWARD. #g"  $FILE > temp; mv temp $FILE

##### forward simulation (data) #####
./bin/xmeshfem2D

if [ $isource -eq 1 ] ; then  
    echo "mpirun -np $NPROC_SPECFEM ./bin/xspecfem2D"
fi
mpirun -np $NPROC_SPECFEM ./bin/xspecfem2D

#### mask source 
# Source location
export xs=$(awk -v "line=$isource" 'NR==line { print $1 }' DATA/sources.dat)
export zs=$(awk -v "line=$isource" 'NR==line { print $2 }' DATA/sources.dat) 
mpirun -np $NPROC_SPECFEM ./bin/mask_func.exe $xs $zs DATA/ OUTPUT_FILES/ 

# save
if [ $isource -eq 1 ]; then  ## for size
    mkdir -p $DISK_DIR/misfit_kernel
    cp -r DATA/proc*_x.bin $DISK_DIR/misfit_kernel/
    cp -r DATA/proc*_z.bin $DISK_DIR/misfit_kernel/
    cp -r DATA/proc*_NSPEC_ibool.bin $DISK_DIR/misfit_kernel/
    cp -r DATA/proc*_jacobian.bin $DISK_DIR/misfit_kernel/
fi
