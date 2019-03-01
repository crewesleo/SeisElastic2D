#!/bin/bash

isource=$1
NPROC_SPECFEM=$2
data_type=$3
velocity_dir=$4
SAVE_FORWARD=$5
WORKING_DIR=$6
DISK_DIR=$7

if [ $isource -eq 1 ] ; then
    echo "SPECFEM3D Adjoint Modeling ..."
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

cp DATA/STATIONS_FILTERED DATA/STATIONS_ADJOINT

##### edit 'Par_file' #####
FILE="./DATA/Par_file"
sed -e "s#^SIMULATION_TYPE.*#SIMULATION_TYPE = 3 #g"  $FILE > temp; mv temp $FILE
sed -e "s#^SAVE_FORWARD.*#SAVE_FORWARD = .$SAVE_FORWARD. #g"  $FILE > temp; mv temp $FILE

##### adjoint simulation  #####
# creates and decomposes mesh
# echo "running mesher..."
mpirun -np $NPROC_SPECFEM ./bin/xmeshfem3D

# runs database generation
echo "running database generation..."
mpirun -np $NPROC_SPECFEM ./bin/xgenerate_databases

# runs simulation
echo "running solver..."
mpirun -np $NPROC_SPECFEM ./bin/xspecfem3D
mv  OUTPUT_FILES/output_solver.txt  OUTPUT_FILES/output_adjoint.txt

# save
cp SEM/*    $ISRC_DATA_DIR/
if [ $isource -eq 1 ]; then  # for size
    echo "copy mesh file to $DISK_DIR/misfit_kernel"
    mkdir -p $DISK_DIR/misfit_kernel 
    cp -r ./OUTPUT_FILES/DATABASES_MPI/proc*_external_mesh.bin $DISK_DIR/misfit_kernel/
fi
