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
    print(qCmd)
    #os.system(qCmd)