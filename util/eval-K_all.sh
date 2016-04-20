#! /bin/bash -x

source ./explist.inc

for e in ${explist[@]}
do
   cd ../${e}/run
      echo; echo ${e}
      echo; echo "msg:"
      more run-K.sh.o*      | grep -A 10 "fapp -C"

      fapppx -A -plimit=1 -Ihwm -d prof -o Profile_fapp.txt || exit 1
      echo; echo "profiler:"
      more Profile_fapp.txt | grep -A 19 "Performance monitor"
   cd -
done
