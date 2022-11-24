#!/bin/sh
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh 
setupATLAS
export ALRB_rootVersion=6.20.06-x86_64-centos7-gcc8-opt
lsetup root 
runNumber=$1
dir=/afs/cern.ch/user/c/ckitsaki/giftests_chara/July2022/final/forGIT/src
cd $dir

root -b -q -l "run.C(\"$runNumber\")"
