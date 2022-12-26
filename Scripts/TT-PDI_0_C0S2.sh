#!/bin/sh
#PBS -l walltime=36:00:00
#PBS -l select=1:ncpus=100:mem=60gb:avx=true

module load gaussian/g16-c01-avx
cp $PBS_O_WORKDIR/TT-PDI_0_C0S2.gjf ./
timeout 36h g16 TT-PDI_0_C0S2.gjf 
formchk TT-PDI_0_C0S2.chk TT-PDI_0_C0S2.fchk
cp *.log  $PBS_O_WORKDIR
    cp *.chk  $PBS_O_WORKDIR
    cp *.fchk  $PBS_O_WORKDIR

    