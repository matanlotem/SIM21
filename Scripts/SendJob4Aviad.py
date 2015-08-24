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

queName = 'barkana'
import os
import datetime

matlabCmd = 'disp(\'Hello World\');'
codePath = r'/a/home/cc/tree/taucc/students/physics/aviadco1/CODE/' #Simulation Runs from here
outputPath = r'/a/home/cc/tree/taucc/students/physics/aviadco1/LOGS/' # Matlab output
outputName = 'Test3' # Matlab output file name
otherArgs = '-l pmem=10gb,pvmem=15gb,nodes=compute-0-62:ppn=1'
sendJob(queName,codePath,outputPath,outputName,matlabCmd,otherArgs)