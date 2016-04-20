#! /bin/bash -x

source ./explist.inc

for e in ${explist[@]}
do
   cd ../${e}/run
      qsub run-linux.sh || exit 1
   cd -
done
