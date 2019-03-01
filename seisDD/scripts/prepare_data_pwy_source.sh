#!/bin/bash

# parameters
source ./parameter
data_tag='DATA_obs'
velocity_dir=$1
SAVE_FORWARD=false

# local id (from 0 to $ntasks-1)
if [ $system == 'slurm' ]; then
    iproc=$SLURM_PROCID  
elif [ $system == 'pbs' ]; then
    iproc=$PBS_VNODENUM
fi

# allocate tasks over all sources
# ntasks in parallel and nsrc in total
# nsrc_per_task=$(( $NSRC / $ntasks ))
# take ceiling 
nsrc_per_task_ceiling=$(echo $(echo "$NSRC $ntasks" | awk '{ print $1/$2 }') | awk '{printf("%d\n",$0+=$0<0?0:0.999)}')
ntasks_ceiling=$(echo $(echo "$NSRC $ntasks" | awk '{print $1%$2}')) 
# take floor 
nsrc_per_task_floor=$(echo $(echo "$NSRC $ntasks" | awk '{ print int($1/$2) }'))

# allocate nsrc for each task
if [ $iproc -lt $ntasks_ceiling ]; then
    nsrc_this_task=$nsrc_per_task_ceiling
    isource_start=$(echo $(echo "$iproc $nsrc_per_task_ceiling" | awk '{ print $1*$2 }'))
else
    nsrc_this_task=$nsrc_per_task_floor
    isource_start=$(echo $(echo "$iproc $nsrc_per_task_floor $ntasks_ceiling $nsrc_per_task_ceiling" | awk '{ print ($1-$3)*$2+$3*$4 }'))
fi

if [ $iproc -eq  0 ]; then
    echo "allocate $NSRC sources over $ntasks tasks"
    echo "iproc 0-$(($ntasks_ceiling-1)): nsrc_per_task=$nsrc_per_task_ceiling"
    echo "iproc $ntasks_ceiling-$(($ntasks-1)): nsrc_per_task=$nsrc_per_task_floor"
fi
echo "iproc = $iproc, isource_start = $isource_start, nsrc_this_task = ${nsrc_this_task}, source: $(($isource_start+1))-$(($isource_start+$nsrc_this_task)) "

# source for this task
for ((isrc_this_task=1; isrc_this_task<=${nsrc_this_task}; isrc_this_task++));
do 
    let isource=`echo $isource_start, $isrc_this_task |awk '{print $1 + $2 }'` 
    if $DISPLAY_DETAILS; then
        echo "iproc = $iproc, isource = $isource"
    fi

    STARTTIME=$(date +%s)
    if  $ExistDATA && [ -d "$DATA_DIR" ]; then
        sh $SCRIPTS_DIR/copy_data_pwy_source.sh $isource $data_tag $data_list $WORKING_DIR $DISK_DIR $DATA_DIR
    else
        sh $SCRIPTS_DIR/Forward_specfem2D_pwy_source.sh $isource $NPROC_SPECFEM $data_tag $data_list \
            $velocity_dir $SAVE_FORWARD $WORKING_DIR $DISK_DIR $DATA_DIR $job 2>./job_info/error_Forward_simulation
    fi
 if [ $isource -eq 1 ] ; then
     ENDTIME=$(date +%s)
     Ttaken=$(($ENDTIME - $STARTTIME))
     echo "Data preparation took $Ttaken seconds"
 fi
done
