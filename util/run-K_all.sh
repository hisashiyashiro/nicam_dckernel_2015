#! /bin/bash -x

source ./explist.inc

for e in ${explist[@]}
do
   cd ../${e}/run
      pjsub run-K.sh || exit 1
   cd -
done
