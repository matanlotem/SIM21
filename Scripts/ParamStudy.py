# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 10:39:50 2015

@author: Matan
"""
import os
import datetime

def sendJob(queName,codePath,outputPath,outputName,matlabCmd,otherArgs=''):
    dt=str(datetime.datetime.now()).replace(' ','_').replace(':','-')[:19]
    outputFile=outputPath + outputName + '_' + dt
    qCmd = 'qsub -q ' + queName + ' -N ' + outputName[:15] + ' -e ' + outputPath + ' -o ' + outputPath + ' ' + otherArgs +\
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

def getCases(dataPath):
    cubeNum = 9
    f = file(dataPath)
    data=[line.split('\t') for line in f.read().splitlines()]
    f.close()
    paramCases = {}
    for line in data[1:]:
        paramCases[int(line[0])] = [cubeNum]+line[1:13]
    return paramCases

def getID(paramCases,caseNum):
    return '_'.join([str(val) for val in paramCases[caseNum][:12]])
    
def getMatlab(paramCases,caseNum):
    return 'RunBackgroundsParam(' + ','.join([str(val) for val in paramCases[caseNum]]) + ');'

def sendCasesJobs(paramCases,caseList):
    otherArgs = '-l pmem=10gb,pvmem=15gb,nodes=compute-0-66:ppn=1'
    for caseNum in caseList:
        outputName = 'ParamStudy_' + str(caseNum) + '_' + getID(paramCases,caseNum)
        matlabCmd = getMatlab(paramCases,caseNum)
        sendJob(queName,codePath,outputPath,outputName,matlabCmd,otherArgs)

def removeData(paramCases,caseList,dataFolders):
    """
    Delete simulation output files by caseID
    """
    for dataFolder in dataFolders:
        cmd = ''
        for caseNum in caseList:
            cmd += 'rm ' + dataFolder[0] + '*' + getID(paramCases,caseNum) + '*\n'
        runOnNodes(dataFolder[1],cmd)
        
def runOnNodes(nodeList,cmd):
    """
    Run command on several nodes.
    If nodeList is empty -> run command localy
    """

    if nodeList == []:
        print 'running localy'
        os.system(cmd)
    else:
        for node in nodeList:
            runOnNode(node,cmd)


def runOnNode(node,cmd):
    """
    Run command on remote node
    """
    
    print 'running on ' + node
    bash =  'ssh ' + node + ' << EOF \n' + cmd + '\nEOF'
    os.system(bash)

queName = 'barkana'
codePath = r'/a/home/cc/tree/taucc/students/physics/matanlotem/Work/CODE/Clean/'
outputPath = r'/a/home/cc/tree/taucc/students/physics/matanlotem/Work/Logs/'
paramDataPath = r'/a/home/cc/tree/taucc/students/physics/matanlotem/Work/ParamStudy/ParamStudy.txt'
#paramDataPath = r'C:\Users\Matan\Work\ParamStudy\ParamStudy.txt'

dataFolders = [['/scratch300/matanlotem/Data/',[]],
               ['/scratch/matanlotem/Data/',['compute-0-66']]]
#               ['/scratch/matanlotem/Data/',['compute-0-62','compute-0-64','compute-0-65','compute-0-66']]]


paramCases = getCases(paramDataPath)
#removeData(paramCases,range(1,24),dataFolders)
sendCasesJobs(paramCases,range(59,65))
