#!/bin/bash
queName="barkana"
codePath="/a/home/cc/tree/taucc/students/physics/matanlotem/Work/CODE/Clean/"
outputPath="/a/home/cc/tree/taucc/students/physics/matanlotem/Work/CODE/Logs/"
outputName="QsubRun"
matlabCmd="RunBackgroundsParam(0,0.05,1,16.5,1,1,0.075,0,0,2,1,0,20);"

./SendJob.sh -q "${queName}" -c "${codePath}" -p "${outputPath}" -o "${outputName}" -m "${matlabCmd}"
