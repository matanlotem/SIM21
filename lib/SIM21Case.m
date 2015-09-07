classdef SIM21Case < matlab.mixin.Copyable
    properties
        name;
        pathExt;
        dataPath;
        tmpDataPath;
        outputPath;
        ncube
        fstar;
        vbc;
        vc;
        fx;
        sed;
        tau;
        feedback;
        delayParam;
        pop;
        fsfunc;
        phVersion;
        zeta;
        ID;
    end


    methods
        function obj = SIM21Case(name,pathExt,outputPath,ncube,fstar,vbc,vc,fx,sed,tau,feedback,delayParam,pop,fsfunc,phVersion,zeta)
            obj.name = name;

            obj.setPathExt(pathExt);
            obj.outputPath = '';
            obj.setOutputPath(outputPath);
            
            obj.ncube = ncube;
            obj.fstar = fstar;
            obj.vbc = vbc;
            obj.vc = vc;
            obj.fx = fx;
            obj.sed = sed;
            obj.tau = tau;
            obj.feedback = feedback;
            obj.delayParam = delayParam;
            obj.pop = pop;
            obj.fsfunc = fsfunc;
            obj.phVersion = phVersion;
            obj.zeta = zeta;
            obj.ID = obj.getID();
        end


        function runID = getID(obj)
            runID = ['_',num2str(obj.ncube),'_',num2str(obj.fstar),'_',num2str(obj.vbc),'_',num2str(obj.vc),'_',num2str(obj.fx),...
                     '_',num2str(obj.sed),'_',num2str(obj.tau),'_',num2str(obj.feedback),'_',num2str(obj.delayParam),'_',num2str(obj.pop),...
                     '_',num2str(obj.fsfunc),'_',num2str(obj.phVersion)];
        end


        function setPathExt(obj,pathExt)
            if ~ isempty(pathExt)
                if pathExt(end) ~= '/'
                    pathExt = [pathExt,'/'];
                end
            end
            obj.pathExt = pathExt;
            obj.dataPath = [SIM21Utils.paths.data,pathExt];
            obj.tmpDataPath = [SIM21Utils.paths.tmpData,pathExt];
        end


        function setOutputPath(obj,newOutputPath,varargin)
            if newOutputPath(end) ~= '/'
                newOutputPath = [newOutputPath,'/'];
            end
            if exist(newOutputPath) ~= 7
                if nargin == 3 && varargin{end}

                    mkdir(newOutputPath);
                    obj.outputPath = newOutputPath;
                else
                    disp (['path doesn''t exist: ',newOutputPath]);
                end
            else
                obj.outputPath = newOutputPath;
            end
        end


        function dataMat = getData(obj,dataType,varargin)
            dataType = SIM21Utils.getDataType(dataType);
            if nargin==3 % If z-range specifed
                zs = varargin{1};
            else % Else - take full z-range
                zs = dataType.z;
            end

            % Create empty matrix
            dataMat = NaN([length(zs),SIM21Utils.cubeSize]);

            % LoadmMatrix for each z if exists
            for zind = 1:length(zs)
                fileName = SIM21Utils.getDataFileName(obj,dataType,zs(zind));
                if exist(fileName) == 2
                    dataMat(zind,:,:,:) = importdata(fileName);
                end
                if mod(zind,10) == 0
                    disp(['    ',num2str(zind),' / ',num2str(length(zs))]);
                end
            end
            dataMat = squeeze(dataMat);
        end


        function runSimulation(obj)
            disp('==Running Full Simulation==');
            curPath = pwd;
            cd(SIM21Utils.paths.code);
            RunBackgroundsParam2(obj.pathExt,obj.ncube,obj.fstar,obj.vbc,obj.vc,obj.fx,obj.sed,obj.tau,obj.feedback,obj.delayParam,obj.pop,obj.fsfunc,obj.phVersion,obj.zeta);
            cd(curPath);
        end


        function runZ(obj,z)
            tic;
            disp(['==Running z=',num2str(z),'==']);
            curPath = pwd;
            cd(SIM21Utils.paths.code);

            global pathname_Data1
            global pathname_Data2
            global pathname_DataBackgrounds
            global delta_cube
            global ID

            pathname_Data1 = obj.tmpDataPath;
            pathname_Data2 = obj.dataPath;
            pathname_DataBackgrounds = SIM21Utils.paths.dataBackgrounds;
            ID = obj.ID;
            delta_cube=importdata(strcat(pathname_DataBackgrounds,'my',num2str(obj.ncube),'_v.dat'));

            BackgroundsParamII(z,obj.ncube,obj.fstar,obj.vbc,obj.vc,obj.fx,obj.sed,obj.tau,obj.zeta,obj.feedback,obj.delayParam,obj.pop,obj.fsfunc,~~obj.phVersion,obj.phVersion);
            
            cd(curPath);
            toc;
        end
    end
end
