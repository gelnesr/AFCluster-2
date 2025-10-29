mkdir -p $SCRATCH/containers $SCRATCH/cache

module load apptainer 
apptainer build $SCRATCH/containers/pytorch_24.12.sif \
  docker://pytorch/pytorch:2.4.1-cuda12.4-cudnn9-runtime

export IMG="$SCRATCH/containers/pytorch_24.12.sif"
export CACHE="$SCRATCH/cache"