classdef findZeta < handle
    properties
        p;
        zetaCases;
    end

    methods
        function obj =findZeta(obj)
            obj.initFindZeta();
            %obj.runJobs(1,[19,20,21]);
            %obj.runJobs(2,[16,21,31]);
            %obj.runJobs(3,[31,41,51]);
            %for caseNum = [10,11,18,19]
            %for caseNum = [1,2,3]
            %    for zeta = obj.getZetasByCase(caseNum)
            %        zetaCase = obj.getZetaCase(caseNum,zeta);
            %        %obj.runJobs(obj.getZetaCase(caseNum,zeta));
            %        obj.getZetaResults(zetaCase);
            %    end
            %end
        end


        function initFindZeta(obj)
            obj.p = paramStudy();

            zetaDataPath = 'findZeta.xlsx';
            rawData = readtable(zetaDataPath);
            for zetaInd = 1:height(rawData)
                zetaCase = obj.initZetaCase(rawData.CASE(zetaInd),rawData.ZETA(zetaInd),rawData.PTAU(zetaInd));
                if zetaInd == 1
                    obj.zetaCases = zetaCase;
                else
                    obj.zetaCases(zetaInd) = zetaCase;
                end
            end
        end


        function zetaCase = initZetaCase(obj,caseNum,zeta,ptau)
            % set paths
            zetaCase.c = copy(obj.p.paramCases(caseNum).c);
            zetaCase.runName = ['findZeta_',num2str(caseNum),'_',num2str(zeta)];
            zetaCase.c.dataPath = ['/scratch300/matanlotem/Data/',zetaCase.runName,'/'];
            zetaCase.c.tmpDataPath = ['/scratch300/matanlotem/TmpData/',zetaCase.runName,'/'];
            zetaCase.c.setOutputPath(['/scratch300/matanlotem/ParamStudy/',zetaCase.runName,'/'],1);

            % set properties
            zetaCase.caseNum = caseNum;
            zetaCase.zeta = zeta;
            zetaCase.c.zeta = zeta;
            zetaCase.ptau = ptau;
            zetaCase.isrun = exist(SIM21Utils.getDataFileName(zetaCase.c,'xHI',6))==2;
            zetaCase = obj.getZetaResults(zetaCase);
        end


        function zetaNum = getZetaNum(obj,caseNum,zeta)
            zetaNum = find(([obj.zetaCases.caseNum] == caseNum).*([obj.zetaCases.zeta] == zeta));
        end


        function zetas = getZetasByCase(obj,caseNum)
            zetas = [obj.zetaCases(find([obj.zetaCases.caseNum] == caseNum)).zeta];
        end
        function zetaCase = getZetaCase(obj,caseNum,zeta)
            zetaCase = obj.zetaCases(obj.getZetaNum(caseNum,zeta));
        end
        function setZetaCase(obj,zetaCase)
            obj.zetaCases(obj.getZetaNum(zetaCase.caseNum,zetaCase.zeta)) = zetaCase;
        end


        function zetaCase = getZetaResults(obj,zetaCase)
            if zetaCase.isrun
                zetaCase.tau = SIM21Analysis.checkTau(zetaCase.c);
                xHIData = SIM21Analysis.getZData(zetaCase.c,'xHI');
                zetaCase.finalXHI = xHIData(2,1);
            else
                zetaCase.tau = NaN;
                zetaCase.finalXHI = NaN;
            end
        end
        function setZetaResults(obj,zetaCase)
            obj.setZetaCase(getZetaResults(zetaCase));
        end


        function res = resultsTable(obj)
            outputMat = [[obj.zetaCases.caseNum]',[obj.zetaCases.zeta]',[obj.zetaCases.finalXHI]',[obj.zetaCases.tau]',[obj.zetaCases.ptau]'];
            res = cell2table(num2cell(outputMat),'variableNames',{'CASE','ZETA','XHI','TAU','PTAU'});
        end
        function saveResultsTable(obj,fileName)
            writetable(obj.resultsTable,fileName);
        end


        function runJobs(obj,zetaCases)
            for zetaCase = zetaCases;
                SIM21Utils.runSimulation(zetaCase.c,zetaCase.runName);
            end
        end
        function runUnrunJobs(obj,caseNums)
            for caseNum = caseNums
                for zeta = obj.getZetasByCase(caseNum)
                    zetaCase = obj.getZetaCase(caseNum,zeta);
                    if ~zetaCase.isrun
                        [caseNum,zeta]
                        SIM21Utils.runSimulation(zetaCase.c,zetaCase.runName);
                    end
                end
            end
        end
        function runJob(obj,caseNum,zeta)
            zetaCase = obj.getZetaCase(caseNum,zeta);
            SIM21Utils.runSimulation(zetaCase.c,zetaCase.runName);
        end
    end
end