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
script=/home/users/akatil/P-ONE/git/PONE_NuTau/BiGauss/likelihoodfit.py
echo "Will use script: " $script

$i3env python $script

date
endsecond=$(date +%s)
echo "End second: " $endsecond
echo "This job took : "`expr $endsecond - $startsecond`" s"
