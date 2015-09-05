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

    properties(Constant)
        baseDataPath = '/scratch300/matanlotem/Data';
        baseTmpDataPath = '/scratch/matanlotem/Data';
        codePath = '/a/home/cc/tree/taucc/students/physics/matanlotem/Work/SIM21/Clean_Fixed/';
        dataBackgrounsPath = '/scratch300/matanlotem/DataBackgrounds_withPlanck/';
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
            obj.dataPath = [obj.baseDataPath,pathExt];
            obj.tmpDataPath = [obj.baseTmpDataPath,pathExt];
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


        function runSimulation(obj)
            disp('==Running Full Simulation==');
            curPath = pwd;
            cd(obj.codePath);
            RunBackgroundsParam2(obj.pathExt,obj.ncube,obj.fstar,obj.vbc,obj.vc,obj.fx,obj.sed,obj.tau,obj.feedback,obj.delayParam,obj.pop,obj.fsfunc,obj.phVersion,obj.zeta);
            cd(curPath);
        end


        function runZ(obj,z)
            tic;
            disp(['==Running z=',num2str(z),'==']);
            curPath = pwd;
            cd(obj.codePath);

            global pathname_Data1
            global pathname_Data2
            global pathname_DataBackgrounds
            global delta_cube
            global ID

            pathname_Data1 = obj.tmpDataPath;
            pathname_Data2 = obj.dataPath;
            pathname_DataBackgrounds = obj.dataBackgrounsPath;
            ID = obj.ID;
            delta_cube=importdata(strcat(pathname_DataBackgrounds,'my',num2str(obj.ncube),'_v.dat'))

            BackgroundsParamII(z,obj.ncube,obj.fstar,obj.vbc,obj.vc,obj.fx,obj.sed,obj.tau,obj.zeta,obj.feedback,obj.delayParam,obj.pop,obj.fsfunc,~~obj.phVersion,obj.phVersion);
            
            cd(curPath);
            toc;
        end
    end
end
