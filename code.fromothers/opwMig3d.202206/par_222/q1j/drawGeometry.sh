#!/bin/bash

gather=$1.su
ascii_geo=$1.txt
bin_geo=$1.bin

sugethw <$gather key=sx,sy,gx,gy|\
    sed 's/sx=//g'| sed 's/sy=//g' | sed 's/gx=//g' |  sed 's/gy=//g'|\
    sed '/^$/d'| sed 's/  *//g'\
    >$ascii_geo

a2b n1=4 <$ascii_geo >$bin_geo

ntr=`wc -l $ascii_geo | awk '{print $1}'`
echo "ntr=$ntr"

#ntr=40
#ntr=4000
xgraph <$bin_geo \
    title="$1, ntr=$ntr"\
    pairs=1.a nplot=$ntr n=2\
    width=800  height=800\
    marksize=6 mark=7 &


