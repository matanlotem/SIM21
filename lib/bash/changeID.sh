#!/bin/bash
args=("$@")
directory=${args[0]}
oldID=${args[1]}
newID=${args[2]}

for file in $( ls ${directory}/*${oldID}.mat ); do
    mv ${file} ${file:0:${#file}-${#oldID}-4}${newID}.mat
done
