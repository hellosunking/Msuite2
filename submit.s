#!/bin/bash
#SBATCH --job-name=build.index
#SBATCH --partition=cpuPartition
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --error=%j.err
#SBATCH --output=%j.out

#module load XXX
cd $SLURM_SUBMIT_DIR

build.index/build.index.sh /lustre/home/ksun/Genomes/hg38/hg38p14.fa.gz /lustre/home/ksun/Genomes/hg38/ncbiRefSeqCurated.20220228.gtf hg38p14 32

#sh wk.sh
#sleep 10d

## submit job
## sbatch XXX.slurm

## check job
## squeue

## delete job
## scancel XXX

