#!/bin/bash
date
startsecond=$(date +%s)

echo "Starting the clsim job"
echo "Argument line : " $@

echo "Starting cvmfs "
eval `/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/setup.sh`
echo "cvfms ran"

i3env=/home/users/akatil/software/V06-01-02/build/env-shell.sh
echo "Will use i3 environment: " ${i3env}
script=/home/users/akatil/P-ONE/git/PONE_NuTau/NuTau_Simulation/step_3_makeHits_pone.py
echo "Will use script: " $script

RUNNUM=$1

GCDFILE=/home/users/akatil/P-ONE/GCD_files/PONE_Phase1.i3.gz
INNAME=step_2_${RUNNUM}_medium_water_NuTau_NuE.i3.gz
INDIR=/data/p-one/akatil/step_2_medium_water/NuTau_NuE_100E_300R_1300H
OUTNAME=step_3_${RUNNUM}_medium_water_custom_mDOM.i3.gz
OUTDIR=/data/p-one/akatil/step_3_medium_water/Custom/NuTau_NuE_100E_300R_1300H

echo "RUNNUM : "$RUNNUM
echo "OUTPUT FILE NAME : "$OUTNAME
echo "OUTPUT FILE DIR  : "$OUTDIR
echo "INPUT FILE NAME : "$INNAME
echo "INPUT FILE DIR  : "$INDIR
echo "GCD FILE : "$GCDFILE

$i3env python $script -n ${RUNNUM} -g ${GCDFILE} -i ${INDIR}/${INNAME} -o ${OUTDIR}/${OUTNAME}

date
endsecond=$(date +%s)
echo "End second: " $endsecond
echo "This job took : "`expr $endsecond - $startsecond`" s"
