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
    qCmd = 'qsub -q ' + queName + ' -N ' + outputName[:15] + ' -e ' + outputPath + ' -o ' + outputPath +\
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
    print(qCmd)
    #os.system(qCmd)

def getCases(dataPath):
    qubeNum = 9
    f = file(dataPath)
    data=[line.split('\t') for line in f.read().splitlines()]
    f.close()
    paramCases = {}
    for line in data[1:]:
        paramCases[int(line[0])] = [qubeNum]+line[1:13]
    return paramCases

def getID(paramCases,caseNum):
    return '_'.join([str(val) for val in paramCases[caseNum][:12]])
    
def getMatlab(paramCases,caseNum):
    return 'RunBackgroundsParam(' + ','.join([str(val) for val in paramCases[caseNum]]) + ');'

def SendCasesJobs(paramCases,caseList):
    for caseNum in caseList:
        outputName = 'ParamStudy_' + str(caseNum) + '_' + getID(paramCases,caseNum)
        matlabCmd = getMatlab(paramCases,caseNum)
        SendJob(queName,codePath,outputPath,outputName,matlabCmd)


queName = 'barkana'
codePath = r'/a/home/cc/tree/taucc/students/physics/matanlotem/Work/CODE/Clean/'
outputPath = r'/a/home/cc/tree/taucc/students/physics/matanlotem/Work/Logs/'

paramDataPath = r'/a/home/cc/tree/taucc/students/physics/matanlotem/Work/ParamStudy/ParamStudy.txt'
paramDataPath = r'C:\Users\Matan\Work\ParamStudy\ParamStudy.txt'

paramCases = getCases(paramDataPath)
SendCasesJobs(paramCases,range(1,9))