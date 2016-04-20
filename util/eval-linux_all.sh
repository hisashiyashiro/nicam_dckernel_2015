#! /bin/bash -x

source ./explist.inc

for e in ${explist[@]}
do
   cd ../${e}/run
      echo; echo ${e}
      cat OUT.log
   cd -
done
