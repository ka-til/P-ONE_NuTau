#!/bin/bash

date
startsecond=$(date +%s)

echo "Starting the NuGen job"
echo "Argument line : " $@

echo "Starting cvmfs "
eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.0.1/setup.sh`
echo "cvfms ran"

i3env=/cvmfs/icecube.opensciencegrid.org/py3-v4.0.1/RHEL_7_x86_64/metaprojects/combo/stable/env-shell.sh
echo "Will use i3 environment: " ${i3env}
script=/home/users/akatil/P-ONE/git/PONE_NuTau/BiGauss/runLikelihoodfit.py
echo "Will use script: " $script

RUNNUM=$1
CMIN=$2
CMAX=$3

INNAME=step_4_${RUNNUM}_medium_water_custom_mDOM_noise.i3.gZ
INDIR=/data/p-one/akatil/step_4_medium_water/NuTau_NuE_20Events
OUTNAME=likelihood_${RUNNUM}
OUTDIR=/data/p-one/akatil/analysis/25_30

$i3env python $script -cmin ${CMIN} -cmax ${CMAX}

date
endsecond=$(date +%s)
echo "End second: " $endsecond
echo "This job took : "`expr $endsecond - $startsecond`" s"
