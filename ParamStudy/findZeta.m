classdef findZeta < handle
    properties
        p;
        zetaCases;
    end

    methods
        function obj =findZeta(obj)
            obj.initFindZeta();
        end


        function initFindZeta(obj)
            obj.p = paramStudy();

            zetaDataPath = 'findZeta.xlsx';
            rawData = readtable(zetaDataPath);
            obj.zetaCases = [];
            for zetaInd = 1:height(rawData)
                obj.zetaCases = [obj.zetaCases, obj.initZetaCase(rawData.CASE(zetaInd),rawData.ZETA(zetaInd),rawData.PTAU(zetaInd),rawData.XHI(zetaInd),rawData.TAU(zetaInd),rawData.SELECTED(zetaInd))];
            end
        end


        function zetaCase = initZetaCase(obj,caseNum,zeta,ptau,finalxHI,tau,selected)
            % set paths
            zetaCase.c = copy(obj.p.paramCases(caseNum).c);
            zetaCase.runName = ['findZeta_',num2str(caseNum),'_',num2str(zeta)];
            zetaCase.c.dataPath = ['/scratch300/matanlotem/Data/',zetaCase.runName,'/'];
            zetaCase.c.tmpDataPath = ['/scratch300/matanlotem/TmpData/',zetaCase.runName,'/'];
            zetaCase.c.setOutputPath(['/scratch300/matanlotem/GRID/',zetaCase.runName,'/']); %%%%%%%%%%%%%

            % set properties
            zetaCase.caseNum = caseNum;
            zetaCase.zeta = zeta;
            zetaCase.c.zeta = zeta;
            zetaCase.ptau = ptau;
            zetaCase.c.ID = zetaCase.c.getID();

            calcResults = 1;
            if exist('finalxHI','var') && exist('tau','var')
                if (~isnan(finalxHI)) && (~isnan(tau))
                    zetaCase.isrun = 1;
                    zetaCase.tau = tau;
                    zetaCase.finalXHI = finalxHI;
                    calcResults = 0;
                end
            end
            if calcResults
                zetaCase.isrun = exist(SIM21Utils.getDataFileName(zetaCase.c,'xHI',6))==2;
                zetaCase = obj.getZetaResults(zetaCase);
            end

            zetaCase.selected = 0;
            if exist('selected','var')
                if selected == 1
                    zetaCase.selected = 1;
                end
            end
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
                disp([zetaCase.caseNum,zetaCase.zeta]);
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
        function runJob(obj,caseNum,zeta,z)
            zetaCase = obj.getZetaCase(caseNum,zeta);
            if exist('z','var')
                SIM21Utils.runSimulation(zetaCase.c,zetaCase.runName,z);
            else
                SIM21Utils.runSimulation(zetaCase.c,zetaCase.runName);
            end
        end


        function copyZetaCase(obj,zetaCase,c,mv)
            cmnd = 'cp';
            if exist('mv','var')
                if mv
                    cmnd = 'mv';
                end
            end
            system(['mkdir ',c.dataPath]);
            system(['mkdir ',c.tmpDataPath]);
            system([cmnd,' ',zetaCase.c.dataPath,'* ',c.dataPath]);
            system([cmnd,' ',zetaCase.c.tmpDataPath,'* ',c.tmpDataPath]);
        end
    end
end