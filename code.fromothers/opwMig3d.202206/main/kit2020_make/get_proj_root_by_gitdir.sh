#!/bin/bash
curr_path=""
found=0
for i in $( seq 1 5 )
do
    curr_targ=$curr_path".git"
    if [ -d $curr_targ ];
    then
        echo $(dirname ${curr_targ})
        found=1
        break
    else
        curr_path="../"$curr_path
    fi
done

if [[ $found == "0" ]]; then
    echo "failed to found .git ! "
fi


