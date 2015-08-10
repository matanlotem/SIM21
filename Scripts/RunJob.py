# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 10:39:50 2015

@author: Matan
"""
import os
import datetime

def SendJob(queName,codePath,outputPath,outputName,matlabCmd):
    dt=str(datetime.datetime.now()).replace(' ','_').replace(':','-')[:19]
    outputFile=outputPath + outputName + '_' + dt
    qCmd = 'qsub -q ' + queName + ' -N ' + outputName + ' -e ' + outputPath + ' -o ' + outputPath +\
           ' << JOBC\n' +\
           'matlab -nodisplay -nosplash -nodesktop -nojvm > ' + outputFile + ' << MATLAB\n' + \
           '    try\n' +\
           '        system(\'hostname\');\n' +\
           '        cd(\'' + codePath + '\');\n' +\
           '        ' + matlabCmd + '\n' +\
           '    catch e\n' +\
           '        disp(e.getReport());\n' +\
           '    end\n' +\
           '    exit;\n' +\
           'MATLAB\n' +\
           'JOBC'
    #print(qCmd)
    os.system(qCmd)


queName = 'barkana'
codePath = '/a/home/cc/tree/taucc/students/physics/matanlotem/Work/CODE/Clean/'
outputPath = '/a/home/cc/tree/taucc/students/physics/matanlotem/Work/Logs/'
outputName = 'RunJob'
matlabCmd = 'RunBackgroundsParam(0,0.05,1,16.5,1,1,0.075,0,0,2,1,0,20);'
SendJob(queName,codePath,outputPath,outputName,matlabCmd)
