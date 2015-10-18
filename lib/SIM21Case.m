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


        function ID = getID(obj)
            ID = ['_',num2str(obj.ncube),'_',num2str(obj.fstar),'_',num2str(obj.vbc),'_',num2str(obj.vc),'_',num2str(obj.fx),...
                  '_',num2str(obj.sed),'_',num2str(obj.tau),'_',num2str(obj.feedback),'_',num2str(obj.delayParam),'_',num2str(obj.pop),...
                  '_',num2str(obj.fsfunc),'_',num2str(obj.phVersion)];
        end


        function setPathExt(obj,pathExt)
            if ~isempty(pathExt)
                if pathExt(end) ~= '/'
                    pathExt = [pathExt,'/'];
                end
            end
            obj.pathExt = pathExt;
            obj.dataPath = [SIM21Utils.paths.data,pathExt];
            obj.tmpDataPath = [SIM21Utils.paths.tmpData,pathExt];
        end


        function setOutputPath(obj,newOutputPath)
            if ~isempty(newOutputPath)
                if newOutputPath(end) ~= '/'
                    newOutputPath = [newOutputPath,'/'];
                end
            end
            if exist(newOutputPath) ~= 7
                mkdir(newOutputPath);
            end
            obj.outputPath = newOutputPath;
        end


        function dataMat = getData(obj,dataType,zs)
            dataType = SIM21Utils.getDataType(dataType);
            if ~ exist('zs','var')
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
            SIM21Utils.runZ(obj,'');
        end


        function runZ(obj,z)
            SIM21Utils.runZ(obj,'',z);
        end
    end
end
