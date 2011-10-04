#!/bin/bash

how_many=2000
for i in {1..40}
    do
        echo $i
        vert=$(($i*$how_many))
        curr_name=`printf "data/sint_distribution/f%03d.png" $i`
        ./main $vert
        python  data/plot.py $curr_name
    done

