classdef SIM21Utils
    properties(Constant)
        % Important Paths
        workPath = '/a/home/cc/tree/taucc/students/physics/matanlotem/Work/';
        paths = struct('work',SIM21Utils.workPath,'code',[SIM21Utils.workPath,'SIM21/Clean_Fixed/'],...
                       'lib',[SIM21Utils.workPath,'SIM21/lib/'],'matrices',[SIM21Utils.workPath,'SIM21/lib/Matrices/'],...
                       'dataBackgrounds','/scratch300/matanlotem/DataBackgrounds_withPlanck/',...
                       'data','/scratch300/matanlotem/Data/','tmpData','/scratch300/matanlotem/TmpData/');

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
                      'nodes','compute-0-66',...
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
            ID = ['_',num2str(MyCube),'_',num2str(MyStar),'_',num2str(MyVBC),'_',num2str(MyVc),...
                  '_',num2str(MyFX),'_',num2str(MySED),'_',num2str(MyTau),'_',num2str(MyFeed),...
                  '_',num2str(DelayParam),'_',num2str(MyPop),'_',num2str(FSfunc),'_',num2str(photoheatingVersion)];
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


        function runSimulation(c,jobName,z)
            runParameters = {c.ncube,c.fstar,c.vbc,c.vc,c.fx,c.sed,c.tau,c.feedback,c.delayParam,c.pop,c.fsfunc,c.phVersion,c.zeta};
            if exist('z','var')
                runParameters = [runParameters,z];
            end

            % Run from console
            if isempty(jobName)
                disp('==Running Full Simulation==');

                curPath = pwd;
                cd(SIM21Utils.paths.code);
                try
                    RunBackgroundsParam3(c.dataPath,c.tmpDataPath,runParameters{:});
                catch e
                    disp(e.getReport());
                end
                cd(curPath);

            % Run as a separated job             
            else
                disp('==Run Simulation Job==');

                resources = 'pmem=11gb,pvmem=17gb';
                if ~ isequal(SIM21Utils.jobs.nodes,'all')
                    resources = [resources,',nodes=',SIM21Utils.jobs.nodes];
                end

                simMatlab = ['RunBackgroundsParam3(''',c.dataPath,''',''',c.tmpDataPath,''',',strjoin(cellfun(@num2str,runParameters,'uniformOutput',false),','),');'];
                SIM21Utils.sendJob(SIM21Utils.paths.code,jobName,simMatlab,resources);
            end
        end


        function runZ(c,z,jobName,comp)
            runParameters = {z,c.ncube,c.fstar,c.vbc,c.vc,c.fx,c.sed,c.tau,c.zeta,c.feedback,c.delayParam,c.pop,c.fsfunc,~~c.phVersion,c.phVersion};

            % Run from console
            if ~ exist('jobName','var')
                disp(['==Running z=',num2str(z),'==']);
                curPath = pwd;
                tic;
                cd(SIM21Utils.paths.code);
                try
                    global pathname_Data1
                    global pathname_Data2
                    global pathname_DataBackgrounds
                    global delta_cube
                    global vbc_cube
                    global ID

                    pathname_Data1 = c.tmpDataPath;
                    pathname_Data2 = c.dataPath;
                    pathname_DataBackgrounds = SIM21Utils.paths.dataBackgrounds;
                    ID = c.ID;
                    delta_cube=importdata(strcat(pathname_DataBackgrounds,'my',num2str(c.ncube),'_d.dat'));
                    vbc_cube=importdata(strcat(pathname_DataBackgrounds,'my',num2str(c.ncube),'_v.dat'));

                    BackgroundsParamII(runParameters{:});
                catch e
                    disp(e.getReport());
                end
                toc;
            cd(curPath);

            % Run as a separated job             
            else
                disp(['==Run z=',num2str(z),' Job==']);

                resources = 'pmem=11gb,pvmem=17gb';
                if exist('comp','var')
                    resources = [resources,',nodes=',comp];
                elseif ~ isequal(SIM21Utils.jobs.nodes,'all')
                    resources = [resources,',nodes=',SIM21Utils.jobs.nodes];
                end

                simMatlab = ['global pathname_Data1\n',...
                             'global pathname_Data2\n',...
                             'global pathname_DataBackgrounds\n',...
                             'global delta_cube\n',...
                             'global vbc_cube\n',...
                             'global ID\n',...
         
                             'pathname_Data1 = ''',c.tmpDataPath,''';\n',...
                             'pathname_Data2 = ''',c.dataPath,''';\n',...
                             'pathname_DataBackgrounds = ',SIM21Utils.paths.dataBackgrounds,';\n',...
                             'ID = ''',c.ID,''';\n',...
                             'delta_cube=importdata(strcat(pathname_DataBackgrounds,''my'',',num2str(c.ncube),',''_d.dat''));\n',...
                             'vbc_cube=importdata(strcat(pathname_DataBackgrounds,''my'',',num2str(c.ncube),',''_v.dat''));\n',...
         
                             'BackgroundsParamII(',strjoin(cellfun(@num2str,runParameters,'uniformOutput',false),','),');'];
                SIM21Utils.sendJob(SIM21Utils.paths.code,jobName,simMatlab,resources);
            end
        end


        function flag = isRun(cases)
            flag = logical([]);
            for ind = 1:length(cases)
                flag(ind) = exist(SIM21Utils.getDataFileName(cases(ind),'xHI',6))==2;
            end
        end
    end
end