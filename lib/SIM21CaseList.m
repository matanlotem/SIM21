classdef SIM21CaseList < handle
	properties
        pathExt;
        outputPath;
        caseList;
    end

    methods
    	function obj = SIM21CaseList(outputPath,pathExt,caseListPath,cubeNum)
    		obj.outputPath = outputPath;
    		obj.pathExt = pathExt;
    		obj.caseList = obj.getCases(caseListPath,cubeNum);
    	end

    	function caseList = getCases(obj,caseListPath,cubeNum)
    		rawData = readtable(caseListPath);
            for i = 1:height(rawData)
                caseNum = rawData.CASE(i);
                caseList(caseNum).caseNum = caseNum;
                caseName = ['Case ',num2str(caseNum)];
                zeta = rawData.Zeta(i);

                if isequal(zeta,'PROBLEM')
                    zeta = NaN;
                    caseList(caseNum).fxHI = NaN;
                    caseList(caseNum).atau = NaN;
                    caseList(caseNum).isgood = false;
                else
                    caseList(caseNum).fxHI = rawData.xHI_z_6_(i);
                    caseList(caseNum).atau = rawData.tau_1(i);
                    caseList(caseNum).isgood = true;
                end
                caseList(caseNum).c = SIM21Case(caseName,[obj.pathExt,num2str(caseNum)],obj.getCaseOutputPath(caseNum),...
                                                cubeNum,rawData.fstar(i),rawData.vbc(i),rawData.vc(i),rawData.fx(i),rawData.sed(i),rawData.tau(i),...
                                                rawData.LWFlag(i),rawData.LWW_s(i),rawData.pop(i),rawData.func(i),rawData.PH(i),zeta);
            end
    	end

    	function caseOutputPath = getCaseOutputPath(obj,caseNum)
            caseOutputPath = [obj.outputPath,'Case_',num2str(caseNum),'/'];
            [status,message,messageid] = mkdir(caseOutputPath);
        end
	end
end
