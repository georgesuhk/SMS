#!/bin/sh
#PBS -l walltime=07:59:00
#PBS -l select=1:ncpus=48:mem=60gb:avx=true

module load gaussian/g16-c01-avx
cp $PBS_O_WORKDIR/TT_1_C1S0.gjf ./
timeout 8h g16 TT_1_C1S0.gjf 
formchk TT_1_C1S0.chk TT_1_C1S0.fchk
cp *.log  $PBS_O_WORKDIR
    cp *.chk  $PBS_O_WORKDIR
    cp *.fchk  $PBS_O_WORKDIR

    