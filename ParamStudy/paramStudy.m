classdef paramStudy < handle
    properties
        dataPath;
        outputPath;
        specialParams;
        paramCases = struct();
    end
    
    
    methods
        function obj = paramStudy()
            addpath('/a/home/cc/tree/taucc/students/physics/matanlotem/Work/SIM21/lib/');
            
            paramDataPath = '/a/home/cc/tree/taucc/students/physics/matanlotem/Work/ParamStudy/ParamStudy.txt';
            cubeNum = 9;
            obj.dataPath = '/scratch300/matanlotem/Data/';
            obj.outputPath = '/scratch300/matanlotem/ParamStudy/';
            
            obj.paramCases = obj.getCases(paramDataPath,cubeNum);
            obj.specialParams = obj.getSpecialParams();
        end
        
        
        function specialParams = getSpecialParams(obj)
            for i=1:length(obj.paramCases)
                i
                if obj.paramCases(i).isgood
                    specialParams(i) = SIM21Analysis.calcSpecialParams(obj.dataPath,obj.outputPath,obj.paramCases(i).ID);
                end
            end
        end


        function paramCases = getCases(obj,paramDataPath,cubeNum)
            rawData = readtable(paramDataPath,'Delimiter','\t');
            paramCases = struct();
            for i = 1:height(rawData)
                paramCases(i).caseNum = rawData.CASE(i);
                paramCases(i).ncube = cubeNum;
                paramCases(i).fstar = rawData.fstar(i);
                paramCases(i).vbc = rawData.vbc(i);
                paramCases(i).vc = rawData.vc(i);
                paramCases(i).fx = rawData.fx(i);
                paramCases(i).sed = rawData.sed(i);
                paramCases(i).tau = rawData.tau(i);
                paramCases(i).feedback = rawData.LWFlag(i);
                paramCases(i).delayParam = rawData.LWW_s(i);
                paramCases(i).pop = rawData.pop(i);
                paramCases(i).fsfunc = rawData.func(i);
                paramCases(i).photoHeating = rawData.PH(i);
                paramCases(i).zeta = char(rawData.Zeta(i));
                if isequal(paramCases(i).zeta,'PROBLEM')
                    paramCases(i).zeta = nan;
                    paramCases(i).isgood = false;
                else
                    paramCases(i).zeta = str2num(paramCases(i).zeta);
                    paramCases(i).isgood = true;
                end
                paramCases(i).ID = SIM21Utils.getID(paramCases(i).ncube,paramCases(i).fstar,paramCases(i).vbc,paramCases(i).vc,...
                                                    paramCases(i).fx,paramCases(i).sed,paramCases(i).tau,paramCases(i).feedback,...
                                                    paramCases(i).delayParam,paramCases(i).pop,paramCases(i).fsfunc,paramCases(i).photoHeating);
            end
        end
    end
end