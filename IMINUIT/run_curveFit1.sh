#!/bin/bash

date
startsecond=$(date +%s)

echo "Starting the job"
echo "Argument line : " $@

echo "Starting cvmfs "
eval `/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/setup.sh`
echo "cvfms ran"

i3env=/home/users/akatil/software/V06-01-02/build/env-shell.sh
echo "Will use i3 environment: " ${i3env}
script=/home/users/akatil/P-ONE/git/PONE_NuTau/IMINUIT/run_curveFit.py
echo "Will use script: " $script

RUNNUM=$1

INNAME=step_5_${RUNNUM}_medium_water_custom_mDOM_recoPulse.i3.gz
INDIR=/data/p-one/akatil/step_5_medium_water/NuTau_NuE_20Events
OUTNAME=step_6_${RUNNUM}_parameters.i3.gz
OUTDIR=/data/p-one/akatil/step_6_analysis/NuTau_NuE_20Events_expGauss_ampRat2

echo "INPUT FILE NAME : "$INNAME
echo "INPUT FILE DIR  : "$INDIR
echo "OUTPUT FILE NAME : "$OUTNAME
echo "OUTPUT FILE DIR  : "$OUTDIR

$i3env python $script  -o ${OUTDIR}/${OUTNAME} -i ${INDIR}/${INNAME}

date
endsecond=$(date +%s)
echo "End second: " $endsecond
echo "This job took : "`expr $endsecond - $startsecond`" s"
