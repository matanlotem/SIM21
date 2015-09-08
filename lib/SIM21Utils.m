classdef SIM21Utils
    properties(Constant)
        % Important Paths
        workPath = '/a/home/cc/tree/taucc/students/physics/matanlotem/Work/';
        paths = struct('work',SIM21Utils.workPath,'code',[SIM21Utils.workPath,'SIM21/Clean_Fixed/'],...
                       'lib',[SIM21Utils.workPath,'SIM21/lib/'],'matrices',[SIM21Utils.workPath,'SIM21/lib/Matrices/'],...
                       'dataBackgrounds','/scratch300/matanlotem/DataBackgrounds_withPlanck/',...
                       'data','/scratch300/matanlotem/Data/','tmpData','/scratch/matanlotem/Data/');

        % Data Matrix Settings
        dataTypes = struct('xHI', struct('z',[6:0.1:15,16:60],'magic','xHI','tmpData',0),...
                           'TK', struct('z',[5:64],'magic','TK','tmpData',0),...
                           'T21cm', struct('z',[6:60],'magic','T21cm','tmpData',0),...
                           'Neut', struct('z',[6:0.1:15,16:60],'magic','Neut','tmpData',0),...
                           'eps', struct('z',[6:60],'magic','eps','tmpData',1),...
                           'xe', struct('z',[5:64],'magic','xe','tmpData',1),...
                           'JLW', struct('z',[6,66],'magic','JLW','tmpData',1),...
                           'Jalpha', struct('z',[6:60],'magic','Jalpha','tmpData',1),...
                           'Lion', struct('z',[6:60],'magic','Lion','tmpData',1));

        cubeSize = [128,128,128];

        jobs = struct('queName','barkana',...
                      'logPath',[SIM21Utils.paths.work,'Logs/']);
    end
    
    methods(Static)
        function M = importMatrix(MName)
            M = importdata(SIM21Utils.getMatrixPath(MName));
        end
        
        
        function MPath = getMatrixPath(MName)
            if length(MName>4)
                if isequal(MName(end-3:end),'.mat')
                    MName = MName(1:end-4)
                end
            end
            MPath=[SIM21Utils.paths.matrices,MName,'.mat'];
        end
        
        
        function ID = getID(MyCube,MyStar,MyVBC,MyVc,MyFX,MySED,MyTau,MyFeed,DelayParam,MyPop,FSfunc,photoheatingVersion)
            ID = ['_' num2str(MyCube)...
                  '_' num2str(MyStar) '_' num2str(MyVBC) '_' num2str(MyVc)...
                  '_' num2str(MyFX) '_' num2str(MySED) '_' num2str(MyTau)...
                  '_' num2str(MyFeed) '_' num2str(DelayParam) '_' num2str(MyPop) '_' num2str(FSfunc) '_' num2str(photoheatingVersion)];
        end


        function dataType = getDataType(dataType)
            if ischar(dataType)
                dataType = getfield(SIM21Utils.dataTypes,dataType);
            end
        end


        function dataFileName = getDataFileName(c,dataType,z)
            % Get raw data matrix file name
            dataType = SIM21Utils.getDataType(dataType);
            if ~ dataType.tmpData
                dataPath = c.dataPath;
            else
                dataPath = c.tmpDataPath;
            end
            dataFileName = [dataPath,dataType.magic,'_',num2str(z),c.ID,'.mat'];
        end


        function compareZ(c1,c2,z)
            disp(['==Checking z=',num2str(z),'==']);            
            SIM21Utils.compareMagic(c1,c2,z,SIM21Utils.dataTypes.xHI);
            SIM21Utils.compareMagic(c1,c2,z-1,SIM21Utils.dataTypes.TK);
            SIM21Utils.compareMagic(c1,c2,z,SIM21Utils.dataTypes.T21cm);
            SIM21Utils.compareMagic(c1,c2,z,SIM21Utils.dataTypes.Neut);
            SIM21Utils.compareMagic(c1,c2,z,SIM21Utils.dataTypes.eps);
            SIM21Utils.compareMagic(c1,c2,z-1,SIM21Utils.dataTypes.xe);
            SIM21Utils.compareMagic(c1,c2,z,SIM21Utils.dataTypes.JLW);
            SIM21Utils.compareMagic(c1,c2,z,SIM21Utils.dataTypes.Jalpha);
            SIM21Utils.compareMagic(c1,c2,z,SIM21Utils.dataTypes.Lion);
        end
        
        
        function compareMagic(c1,c2,z,dataType)
            if max(find(dataType.z==z))
                msg = '';
                a1 = SIM21Utils.getDataFileName(c1,dataType,z);
                a2 = SIM21Utils.getDataFileName(c2,dataType,z);
                if exist(a1, 'file') ~= 2
                    msg = [msg,'No File'];
                elseif exist(a2, 'file') ~= 2
                    msg = [msg,'No tmp File'];
                elseif isequal(importdata(a1),importdata(a2))
                    msg = [msg,'1'];
                else
                    msg = [msg,'0'];
                end
                disp(sprintf([dataType.magic,'\t',msg]));
            end
        end


        function sendJob(codePath,jobName,matlabCmd,resources)
            outputFile = [SIM21Utils.jobs.logPath,jobName,'_',datestr(datetime,'yyyy-mm-dd_HH-MM-SS')];
            qCmd = sprintf(['qsub -q ',SIM21Utils.jobs.queName,' -N ',jobName(1:min(end,15)),' -e /dev/null -o /dev/null -l ',resources,...
                            ' << JOBC\n',...
                            'matlab -nodisplay -nosplash -nodesktop -nojvm > ',outputFile,' << MATLAB\n',...
                            '    try\n',...
                            '        system(''hostname'');\n',...
                            '        cd(''',codePath,''');\n',...
                            '        ',matlabCmd,'\n',...
                            '    catch e\n',...
                            '        disp(e.getReport());\n',...
                            '    end\n',...
                            '    exit;\n',...
                            'MATLAB\n',...
                            'JOBC']);
            [a,hostname] = system('hostname');
            if ~ isequal(hostname(1:5),'power')
                qCmd = sprintf(['ssh power << ROP\n',qCmd,'\nROP\nexit']);
            end
            %disp(qCmd);
            system(qCmd);
        end


        function runSimJob(c,jobName)
            %resources = 'pmem=10gb,pvmem=15gb,nodes=compute-0-66:ppn=1';
            resources = 'pmem=11gb,pvmem=17gb';
            simMatlab = ['RunBackgroundsParam2(''',c.dataPath,''',''',c.tmpDataPath,''',',...
                                               num2str(c.ncube),',',num2str(c.fstar),',',num2str(c.vbc),',',num2str(c.vc),',',...
                                               num2str(c.fx),',',num2str(c.sed),',',num2str(c.tau),',',num2str(c.feedback),',',num2str(c.delayParam),',',...
                                               num2str(c.pop),',',num2str(c.fsfunc),',',num2str(c.phVersion),',',num2str(c.zeta),');'];
            SIM21Utils.sendJob(SIM21Utils.paths.code,jobName,simMatlab,resources)
        end
    end
end