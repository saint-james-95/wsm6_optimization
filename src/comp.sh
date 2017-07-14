#!/bin/bash
./makewsm6 $1 $2 $3 $4 $5 $6 $7 >& make.xeon.out
if [[ $(grep -i error make.xeon.out) ]]; then
	echo "Error in compiling. See make.xeon.out for details."
	exit -1
fi
SRCDIR=$(grep "SRCDIR = " make.xeon.out | sed 's|.*= \(.*\)|\1|')
RUNDIR=$(grep -A1 "Run scripts are in:" make.xeon.out | tail -1)
RUNSCRPT=$(grep -A1 "Run script is" make.xeon.out | tail -1 | sed 's/^\s*//')
echo $RUNSCRPT
cd $RUNDIR
./$RUNSCRPT
exit 0
