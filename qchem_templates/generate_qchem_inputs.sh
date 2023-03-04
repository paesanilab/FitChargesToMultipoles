#!/bin/bash

if [ $# != 1 ]; then
  echo "USAGE: $0 <configs.xyz>"
  exit
fi

BASEDIR=$(dirname "$0")

nat=`head -n 1 $1 | awk '{print $1}'`
nlxf=$(($nat + 2))
nl=`cat $1 | wc -l`
n=$(($nl/$nlxf))

for i in `seq 1 1 $n`; do
  head -n $(($nlxf*$i)) $1 | tail -n $nat > tmp
  cat $BASEDIR/HEAD tmp $BASEDIR/TAIL > input$i
done
