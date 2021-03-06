#! /bin/bash -x
#
# for K computer
#
#PJM --rsc-list "rscgrp=micro"
#PJM --rsc-list "node=1"
#PJM --rsc-list "elapse=00:30:00"
#PJM -j
#PJM -s
#
. /work/system/Env_base
#
export PARALLEL=8
export OMP_NUM_THREADS=8
export XOS_MMM_L_ARENA_FREE=2
#export fu07bf=1

HMDIR=`pwd`/../..

OUTDIR=${HMDIR}/dckernel_vertical_adv_limiter/run
mkdir -p ${OUTDIR}
cd       ${OUTDIR}

ln -svf ${HMDIR}/bin/dckernel_vertical_adv_limiter.exe .
ln -svf ${HMDIR}/dckernel_vertical_adv_limiter/data/snapshot.Vertical_Adv_limiter.001 .

rm -rf ./prof*

fapp -C -Ihwm -Hevent=Statistics -d prof -L 10 ./dckernel_vertical_adv_limiter.exe
