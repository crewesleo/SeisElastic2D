#!/bin/bash

isource=$1
NPROC_SPECFEM=$2
compute_adjoint=$3
data_list=$4
measurement_list=$5
misfit_type_list=$6
WORKING_DIR=$7
DISK_DIR=$8
Wscale=$9
wavelet_path=${10}
VISCOELASTIC=${11}
measurement_attenuation=${12}

if [ $isource -eq 1 ]; then
    echo "adjoint source ..."
    echo "NPROC_SPECFEM=$NPROC_SPECFEM"
    echo "compute_adjoint=$compute_adjoint"
    echo "data_list=$data_list"
    echo "measurement_list=$measurement_list"
    echo "misfit_type_list=$misfit_type_list"
    echo "Wscale=$Wscale"
    echo "wavelet_path=$wavelet_path"
    echo "WORKING_DIR=$WORKING_DIR"
    echo "DISK_DIR=$DISK_DIR"
    echo "VISCOELASTIC=$VISCOELASTIC"
    echo "measurement_attenuation=$measurement_attenuation"
fi

ISRC_WORKING_DIR=$( seq --format="$WORKING_DIR/%06.f/" $(($isource-1)) $(($isource-1)) )

mkdir -p $ISRC_WORKING_DIR
cd $ISRC_WORKING_DIR

if [ $Wscale -gt 0 ]; then
    cp -r $wavelet_path ./
fi

INPUT_DIR=$( seq --format="$DISK_DIR/%06.f/" $(($isource-1)) $(($isource-1)) )
mkdir -p $INPUT_DIR/SEM_ATTENUATION
#mkdir -p $INPUT_DIR/SEM_WD

# adjoint source
mpirun -np $NPROC_SPECFEM ./bin/misfit_adjoint_attenuation.exe $compute_adjoint $data_list $measurement_list $misfit_type_list $INPUT_DIR $VISCOELASTIC $measurement_attenuation

## copy and postprocessing of adjoint source
arr=$(echo $data_list | tr "," "\n")

for x in $arr
do # commented by PWY 14-01-2018
    if [ -f "SU_process/process_adj.sh" ]; then
        sh SU_process/process_adj.sh \
            $INPUT_DIR/SEM_ATTENUATION/U${x}_file_single.su.adj \
            $ISRC_WORKING_DIR/SEM/U${x}_file_single.su.adj            
    else
        cp  $INPUT_DIR/SEM_ATTENUATION/U${x}_file_single.su.adj \
            $ISRC_WORKING_DIR/SEM/U${x}_file_single.su.adj  
    fi
done

