#!/bin/bash

# parameters
source ./parameter
velocity_dir=$1
compute_adjoint=$2

# local id (from 0 to $ntasks-1)
if [ $system == 'slurm' ]; then
    iproc=$SLURM_PROCID  
elif [ $system == 'pbs' ]; then
    iproc=$PBS_VNODENUM
fi

# allocate tasks over all sources
# ntasks in parallel and nsrc in total
# take ceiling 
nsrc_per_task_ceiling=$(echo $(echo "$NSRC $ntasks" | awk '{ print $1/$2 }') | \
    awk '{printf("%d\n",$0+=$0<0?0:0.999)}')
ntasks_ceiling=$(echo $(echo "$NSRC $ntasks" | awk '{print $1%$2}'))
# take floor 
nsrc_per_task_floor=$(echo $(echo "$NSRC $ntasks" | awk '{ print int($1/$2) }'))


# allocate nsrc for each task
if [ $iproc -lt $ntasks_ceiling ]; then
    nsrc_this_task=$nsrc_per_task_ceiling
    isource_start=$(echo $(echo "$iproc $nsrc_per_task_ceiling" | awk '{ print $1*$2 }'))
else
    nsrc_this_task=$nsrc_per_task_floor
    isource_start=$(echo $(echo "$iproc $nsrc_per_task_floor \
        $ntasks_ceiling $nsrc_per_task_ceiling" | awk '{ print ($1-$3)*$2+$3*$4 }'))
fi

if $DISPLAY_DETAILS; then
    echo "iproc = $iproc, isource_start = $isource_start, nsrc_this_task = ${nsrc_this_task}, \
        source: $(($isource_start+1))-$(($isource_start+$nsrc_this_task)) "
fi

# source for this task
for ((isrc_this_task=1; isrc_this_task<=${nsrc_this_task}; isrc_this_task++));
do
    let isource=`echo $isource_start, $isrc_this_task |awk '{print $1 + $2 }'`
    if $DISPLAY_DETAILS; then
        echo "iproc = $iproc, isource = $isource"
    fi

    # STEP one -- forward simulation
    STARTTIME=$(date +%s)
    data_tag='DATA_syn'
    if $compute_adjoint ; then   
        SAVE_FORWARD=true
    else
        SAVE_FORWARD=false
    fi
    sh $SCRIPTS_DIR/Forward_specfem2D_pwy.sh $isource $NPROC_SPECFEM $data_tag $data_list \
        $velocity_dir $SAVE_FORWARD $WORKING_DIR $DISK_DIR $DATA_DIR $job 2>./job_info/error_Forward_simulation
     if [ $isource -eq 1 ] && $compute_adjoint ; then
         ENDTIME=$(date +%s)
         Ttaken=$(($ENDTIME - $STARTTIME))
         echo "Forward simulation took $Ttaken seconds"
     fi

    # STEP two -- adjoint source
    STARTTIME=$(date +%s)
    sh $SCRIPTS_DIR/adjoint_source.sh $isource $NPROC_SPECFEM $compute_adjoint $data_list \
        $measurement_list $misfit_type_list $WORKING_DIR $DISK_DIR $Wscale $wavelet_path 2>./job_info/error_adj_source
     if [ $isource -eq 1 ] && $compute_adjoint ; then
         ENDTIME=$(date +%s)
         Ttaken=$(($ENDTIME - $STARTTIME))
         echo "adjoint source took $Ttaken seconds"
     fi


    # STEP three -- adjoint simulation?
    STARTTIME=$(date +%s)
    if $compute_adjoint; then
        data_tag='SEM'
        SAVE_FORWARD=false
        sh $SCRIPTS_DIR/Adjoint_${solver}.sh $isource $NPROC_SPECFEM $data_tag \
            $velocity_dir $SAVE_FORWARD $WORKING_DIR $DISK_DIR 2>./job_info/error_Adjoint_simulation
    fi
     if [ $isource -eq 1 ] && $compute_adjoint ; then
         ENDTIME=$(date +%s)
         Ttaken=$(($ENDTIME - $STARTTIME))
         echo "Adjoint simulation took $Ttaken seconds"
     fi

done
