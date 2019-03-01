#!/bin/bash
#SBATCH -J K_att
#SBATCH -N 32
#SBATCH -n 32
#SBATCH -o slurm.log
#SBATCH -t 00:30:00
#SBATCH -A w18_trust

#SBATCH --qos=interactive

ulimit -s unlimited

source ./parameter
export user=$(whoami)

if [ $system == 'slurm' ]; then
    # Submit directory
    export SUBMIT_DIR=$SLURM_SUBMIT_DIR
    echo "$SLURM_JOB_NODELIST"  >  ./job_info/NodeList
    echo "$SLURM_JOBID"  >  ./job_info/JobID
elif [ $system == 'pbs' ]; then
    # Submit directory
    export SUBMIT_DIR=$PBS_O_WORKDIR
    echo "$PBS_NODEFILE"  >  ./job_info/NodeList
    echo "$PBS_JOBID"  >  ./job_info/JobID
fi
cd $SUBMIT_DIR
#################### input parameters ###################################################
# directories
export SCRIPTS_DIR="$package_path/scripts" 
export SUBMIT_RESULT="$SUBMIT_DIR/RESULTS/$job/Scale${Wscale}_${measurement_list}_${misfit_type_list}"  # final results
if [ -z "$working_path" ]; then
   export working_path=$SUBMIT_DIR
fi
export WORKING_DIR="$working_path/$Job_title/specfem/"  # directory on local nodes, where specfem runs
export DISK_DIR="$working_path/$Job_title/output/"      # temporary directory for data/model/gradient ...

echo "Submit job << $Job_title >> in : $SUBMIT_DIR  "
echo "Working directory: $WORKING_DIR"
echo "FINAL results in :  $SUBMIT_RESULT"

#########################################################################################
STARTTIME=$(date +%s)
echo "start time is :  $(date +"%T")"

rm -rf $WORKING_DIR
mkdir -p $WORKING_DIR

if $ReStart; then
    echo
    echo "Re-Starting job ..." 
    echo "Clean up result/DISK directories ..."
    rm -rf $SUBMIT_RESULT $DISK_DIR
    mkdir -p $SUBMIT_RESULT $DISK_DIR
else
    echo
    echo "Continue with current job ..."
fi 

echo 
echo "prepare data ..."
velocity_dir=$target_velocity_dir
if [ $system == 'slurm' ]; then
    srun -n 32 -c $NPROC_SPECFEM -l -W 0 $SCRIPTS_DIR/prepare_data_pwy_source.sh $velocity_dir 2> ./job_info/error_target
elif [ $system == 'pbs' ]; then 
    pbsdsh -n $ntasks -c $NPROC_SPECFEM -l -W 0 $SCRIPTS_DIR/prepare_data.sh $velocity_dir 2> ./job_info/error_target
fi

echo
echo "prepare starting model ..."
rm -rf $DISK_DIR/m_current
cp -r $initial_velocity_dir    $DISK_DIR/m_current

echo
echo "********************************************************************************************************"
echo "       Welcome job << $job >> " 
echo "       Scale: '$Wscale'; measurement: '${measurement_list}'; misfit_type: '${misfit_type_list}' " 
echo "********************************************************************************************************"

echo "Forward/Adjoint simulation for current model ...... "
velocity_dir=$DISK_DIR/m_current
compute_adjoint=true
if [ $system == 'slurm' ]; then
    srun -n 32 -c $NPROC_SPECFEM -l -W 0 $SCRIPTS_DIR/Adjoint_pwy_source.sh $velocity_dir $compute_adjoint 2> ./job_info/error_current
elif [ $system == 'pbs' ]; then
    pbsdsh -n $ntasks -c $NPROC_SPECFEM -l -W 0 $SCRIPTS_DIR/Adjoint.sh $velocity_dir $compute_adjoint 2> ./job_info/error_current
fi
echo 
echo "sum event kernel ...... "
mkdir -p $DISK_DIR/misfit_kernel
mpirun -np $NPROC_SPECFEM ./bin/sum_kernel.exe $kernel_list,$precond_list $WORKING_DIR $DISK_DIR 2> ./job_info/error_sum_kernel

if $VISCOELASTIC ; then
    # remove attenuation kernels first
    rm -r $DISK_DIR/misfit_kernel/*Qkappa*
    rm -r $DISK_DIR/misfit_kernel/*Qmu*
    echo "compute attenuation kernels"
    if [ $system == 'slurm' ]; then
        srun -n 32 -c $NPROC_SPECFEM -l -W 0 $SCRIPTS_DIR/Adjoint_pwy_source_attenuation.sh $velocity_dir $compute_adjoint 2> ./job_info/error_current
    elif [ $system == 'pbs' ]; then
        pbsdsh -n $ntasks -c $NPROC_SPECFEM -l -W 0 $SCRIPTS_DIR/Adjoint_attenuation.sh $velocity_dir $compute_adjoint 2> ./job_info/error_current
    fi
    echo "summing attenuation kernels..."
    mpirun -np $NPROC_SPECFEM ./bin/sum_kernel.exe 'Qkappa_kernel,Qmu_kernel' $WORKING_DIR $DISK_DIR 2> ./job_info/error_sum_kernel_attenuation
fi

if $smooth ; then
    echo 
    echo "smooth misfit kernel ... "
    if [ $solver == 'specfem3D' ]; then 
        rm -rf OUTPUT_FILES 
        mkdir OUTPUT_FILES
        mkdir OUTPUT_FILES/DATABASES_MPI
        cp $DISK_DIR/misfit_kernel/proc*external_mesh.bin OUTPUT_FILES/DATABASES_MPI/   
    fi
    mpirun -np $NPROC_SPECFEM ./bin/xsmooth_sem $sigma_x $sigma_z $kernel_list,$precond_list $DISK_DIR/misfit_kernel/ $DISK_DIR/misfit_kernel/ $GPU_MODE 2> ./job_info/error_smooth_kernel
fi

echo
echo "******************finish all for scale $Wscale **************"

cp -r $SUBMIT_DIR/parameter $SUBMIT_RESULT/
cp -r $DISK_DIR/misfit_kernel $SUBMIT_RESULT/
#mkdir -p $SUBMIT_RESULT/WD
#cp -r $DISK_DIR/misfit_kernel $SUBMIT_RESULT/WD
echo
echo " clean up local nodes (wait) ...... "
#rm -rf $working_path/$Job_title
#rm -rf OUTPUT_FILES

ENDTIME=$(date +%s)
Ttaken=$(($ENDTIME - $STARTTIME))
echo
echo "finish time is : $(date +"%T")" 
echo "RUNTIME is :  $(($Ttaken / 3600)) hours ::  $(($(($Ttaken%3600))/60)) minutes  :: $(($Ttaken % 60)) seconds."

echo
echo "******************well done*******************************"

cp -r $SUBMIT_DIR/job_info/output $SUBMIT_RESULT/
