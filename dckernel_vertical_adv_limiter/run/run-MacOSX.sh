#! /bin/bash -x

HMDIR=`pwd`/../..

OUTDIR=${HMDIR}/dckernel_vertical_adv_limiter/run
mkdir -p ${OUTDIR}
cd       ${OUTDIR}

ln -svf ${HMDIR}/bin/dckernel_vertical_adv_limiter.exe .
ln -svf ${HMDIR}/dckernel_vertical_adv_limiter/data/snapshot.Vertical_Adv_limiter.001 .

./dckernel_vertical_adv_limiter.exe
