#!/bin/bash
nodeNums=(62 64 65 66)
for nodeNum in ${nodeNums[@]};do
    ssh compute-0-$nodeNum << EOF
    mkdir /scratch/matanlotem
    mkdir /scratch/matanlotem/Data
    exit
EOF
done
