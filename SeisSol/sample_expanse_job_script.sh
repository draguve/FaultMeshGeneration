#!/bin/bash

# Output and error (also --output, --error):
#SBATCH -o ./%j.%x.out
#SBATCH -e ./%j.%x.err
#SBATCH --job-name="Seissol-double"
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ritwik.patil+slurm@sjsu.edu
# Wall clock limit:
#SBATCH --time=23:00:00
#SBATCH --no-requeue
#SBATCH --partition=compute
#SBATCH --account=ddp454
#SBATCH --mem=0
#SBATCH --constraint="lustre"
# Setup of execution environment
#SBATCH --export=ALL
# Number of nodes and MPI tasks per node:
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --exclusive

# Run the program:
export MP_SINGLE_THREAD=yes
unset KMP_AFFINITY
export OMP_NUM_THREADS=127
export OMP_PLACES="cores(127)"
export XDMFWRITER_ALIGNMENT=8388608
export XDMFWRITER_BLOCK_SIZE=8388608
export SC_CHECKPOINT_ALIGNMENT=8388608
export SEISSOL_CHECKPOINT_ALIGNMENT=8388608
export SEISSOL_CHECKPOINT_DIRECT=1
export ASYNC_MODE=THREAD
export ASYNC_BUFFER_ALIGNMENT=8388608

echo 'num_nodes:' ${SLURM_JOB_NUM_NODES} 'ntasks:' ${SLURM_NTASKS}

# harness script will find the output directory and make it if it does not exist, store into output directory variable
export OUTPUT_DIR=$(python ~/FaultMeshGeneration/SeisSol/harness.py before)
export CONFIG_DIR=$(dirname "$OUTPUT_DIR")
export FOLDER_NAME=$(basename "$CONFIG_DIR")

# start simulation
ulimit -Ss 2097152
srun --mpi=pmix_v3 -n ${SLURM_NTASKS} SeisSol_Release_drome_4_elastic parameters.par

# Postprocess will generate the plots from the energy file and print out useful metrics
export NUMEXPR_NUM_THREADS=16
python ~/FaultMeshGeneration/SeisSol/harness.py after

python ~/FaultMeshGeneration/SeisSol/thenotetaker.py

# If output folder isn't empty then compress it
if [ -d "${OUTPUT_DIR}" ] && [ "$(ls -A "${OUTPUT_DIR}")" ]; then
    set +e
    tar --use-compress-program="pigz" -cvf "../${FOLDER_NAME}.tar.gz" -C "$(dirname "$CONFIG_DIR")" "${FOLDER_NAME}"
    exitcode=$?

    if [ "$exitcode" != "1" ] && [ "$exitcode" != "0" ]; then
        exit $exitcode
    fi
    set -e
fi
