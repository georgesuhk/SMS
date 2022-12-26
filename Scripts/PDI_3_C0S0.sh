#!/bin/sh
#PBS -l walltime=36:00:00
#PBS -l select=1:ncpus=100:mem=60gb:avx=true

module load gaussian/g16-c01-avx
cp $PBS_O_WORKDIR/PDI_3_C0S0.gjf ./
timeout 36h g16 PDI_3_C0S0.gjf 
formchk PDI_3_C0S0.chk PDI_3_C0S0.fchk
cp *.log  $PBS_O_WORKDIR
    cp *.chk  $PBS_O_WORKDIR
    cp *.fchk  $PBS_O_WORKDIR

    