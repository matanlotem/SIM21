classdef paramStudy < handle
    properties
        pathExt;
        outputPath;
        regularCase;
        smallVarCases;
        largeVarCases;
        otherCases;
        workCases;
        specialParams;
        tempParams;
        paramCases;
    end
    
    
    methods
        function obj = paramStudy()
            addpath('/a/home/cc/tree/taucc/students/physics/matanlotem/Work/SIM21/lib/');
            
            paramDataPath = 'paramStudy.xlsx';
            cubeNum = 9;
            obj.pathExt = 'PS_';
            obj.outputPath = '/scratch300/matanlotem/ParamStudy/';
            obj.regularCase = [1];
            obj.smallVarCases = [2:33];
            obj.largeVarCases = [34:65];
            obj.otherCases = [66:69];
            
            obj.paramCases = obj.getCases(paramDataPath,cubeNum);
            obj.workCases = obj.areRun([obj.regularCase,obj.smallVarCases,obj.largeVarCases,obj.otherCases]);
            obj.specialParams = obj.getSpecialParams();
            obj.tempParams = obj.getTempParams();
            %obj.plotSomething3();
        end


        function paramCases = getCases(obj,paramDataPath,cubeNum)
            rawData = readtable(paramDataPath);
            for i = 1:height(rawData)
                caseNum = rawData.CASE(i);
                paramCases(caseNum).caseNum = caseNum;
                caseName = ['Case ',num2str(caseNum)];
                zeta = rawData.Zeta(i);

                if isequal(zeta,'PROBLEM')
                    zeta = NaN;
                    paramCases(caseNum).fxHI = NaN;
                    paramCases(caseNum).atau = NaN;
                    paramCases(caseNum).isgood = false;
                else
                    paramCases(caseNum).fxHI = rawData.xHI_z_6_(i);
                    paramCases(caseNum).atau = rawData.tau_1(i);
                    paramCases(caseNum).isgood = true;
                end
                paramCases(caseNum).c = SIM21Case(caseName,[obj.pathExt,num2str(caseNum)],obj.getCaseOutputPath(caseNum),...
                                                  cubeNum,rawData.fstar(i),rawData.vbc(i),rawData.vc(i),rawData.fx(i),rawData.sed(i),rawData.tau(i),...
                                                  rawData.LWFlag(i),rawData.LWW_s(i),rawData.pop(i),rawData.func(i),rawData.PH(i),zeta);
            end
        end


        function runCaseNums = areRun(obj,caseNums)
            runCaseNums = caseNums(SIM21Utils.isRun([obj.paramCases(caseNums).c]));
        end
        
        
        function specialParams = getSpecialParams(obj)
            for caseNum = obj.workCases
                paramMatName = [obj.paramCases(caseNum).c.outputPath,'specialParams',obj.paramCases(caseNum).c.ID,'.mat'];
                
                % Check if specialParams exists and load
                if exist(paramMatName) == 2
                    specialParams(caseNum)=importdata(paramMatName);
                else
                    disp(caseNum);
                    specialParams(caseNum) = SIM21Analysis.calcSpecialParams(obj.paramCases(caseNum).c);
                end
            end
        end
        
        
        function recalcSpecialParams(obj)
            % Recalculate specialParams
            for caseNum = obj.workCases
                disp(caseNum);
                specialParams(caseNum) = SIM21Analysis.calcSpecialParams(obj.paramCases(caseNum).c);
            end
            obj.specialParams = specialParams;
        end
        

        function tempParams = getTempParams(obj)
            % Calculate tempParams
            for caseNum = obj.workCases
                tempParams(caseNum).minMaxSlope.zDiff = obj.specialParams(caseNum).maxSlope.z - obj.specialParams(caseNum).minSlope.z;
                tempParams(caseNum).xHI75_25.zDiff = obj.specialParams(caseNum).xHI75.z - obj.specialParams(caseNum).xHI25.z;
                tempParams(caseNum).vcFstar = obj.paramCases(caseNum).c.vc * obj.paramCases(caseNum).c.fstar;
                tempParams(caseNum).xHI75zNorm = obj.specialParams(caseNum).xHI75.z / obj.paramCases(caseNum).atau;
            end
        end
        

        function specialParamsTable(obj)
            headers = ['General\t\t'   ,'minT21cm\t\t' ,'maxT21cm\t\t' ,'minSlope\t\t\t'      ,'maxSlope\t\t\t'      ,'minTK\t\t' ,'0 Crossing\t\t' ,'xHI percentage\t\t\t\t'       ,'THT\t\t'   ,'final xHI\t\t'        ,'tau\t\t'              ,'\n',...
                       'Case\t','ID\t' ,'z\t','T\t'    ,'z\t','T\t'    ,'z\t','T\t','slope\t' ,'z\t','T\t','slope\t' ,'z\t','T\t' ,'z\t','NOC\t'    ,'z75\t','z50\t','z25\t','z0\t' ,'z\t','T\t' ,'planned\t','actual\t' ,'planned\t','actual\t' ,'\n'];
            tableStr = headers;
            
            function cellVal = formatCell(cellVal)
                if isfloat(cellVal)
                    cellVal = regexprep(num2str(cellVal),'\s*',',');
                end
            end
            
            for caseNum = obj.workCases
                disp(caseNum);
                xHIData = SIM21Analysis.getZData(obj.paramCases(caseNum).c,SIM21Utils.dataTypes.xHI);
                fxHI = xHIData(2,1);
                xHIData = [];
                atau = SIM21Analysis.checkTau(obj.paramCases(caseNum).c);
                
                rowStr = [formatCell(caseNum)                                     ,'\t',...
                          formatCell(obj.paramCases(caseNum).c.ID)                ,'\t',...
                          formatCell(obj.specialParams(caseNum).minT21cm.z)       ,'\t',...
                          formatCell(obj.specialParams(caseNum).minT21cm.T)       ,'\t',...
                          formatCell(obj.specialParams(caseNum).maxT21cm.z)       ,'\t',...
                          formatCell(obj.specialParams(caseNum).maxT21cm.T)       ,'\t',...
                          formatCell(obj.specialParams(caseNum).minSlope.z)       ,'\t',...
                          formatCell(obj.specialParams(caseNum).minSlope.T)       ,'\t',...
                          formatCell(obj.specialParams(caseNum).minSlope.slope)   ,'\t',...
                          formatCell(obj.specialParams(caseNum).maxSlope.z)       ,'\t',...
                          formatCell(obj.specialParams(caseNum).maxSlope.T)       ,'\t',...
                          formatCell(obj.specialParams(caseNum).maxSlope.slope)   ,'\t',...
                          formatCell(obj.specialParams(caseNum).minTK.z)          ,'\t',...
                          formatCell(obj.specialParams(caseNum).minTK.T)          ,'\t',...
                          formatCell(obj.specialParams(caseNum).xCross.z)         ,'\t',...
                          formatCell(length(obj.specialParams(caseNum).xCross.z)) ,'\t',...
                          formatCell(obj.specialParams(caseNum).xHI75.z)          ,'\t',...
                          formatCell(obj.specialParams(caseNum).xHI50.z)          ,'\t',...
                          formatCell(obj.specialParams(caseNum).xHI25.z)          ,'\t',...
                          formatCell(obj.specialParams(caseNum).xHI0.z)           ,'\t',...
                          formatCell(obj.specialParams(caseNum).THT.z)            ,'\t',...
                          formatCell(obj.specialParams(caseNum).THT.T)            ,'\t',...
                          formatCell(obj.paramCases(caseNum).fxHI)                ,'\t',...
                          formatCell(fxHI)                                        ,'\t',...
                          formatCell(obj.paramCases(caseNum).atau)                ,'\t',...
                          formatCell(atau)                                        ,'\n'];
                tableStr = [tableStr,rowStr];
            end
            
            fileName = [obj.outputPath,'specialParams.xls'];
            fid = fopen(fileName,'w');
            fwrite(fid,sprintf(tableStr));
            fclose(fid);
        end
        
        
        function plotSanityGraphs(obj,caseNums)
            if exist('caseNums','var')
                caseNums = intersect(caseNums,obj.workCases);
            else
                caseNums = obj.workCases;
            end
            for caseNum = caseNums
                disp(caseNum);
                SIM21Analysis.plotSanityGraphs(obj.paramCases(caseNum).c);
            end
        end

        function plotAllZGraphs(obj,caseNums)
            if exist('caseNums','var')
                caseNums = intersect(caseNums,obj.workCases);
            else
                caseNums = obj.workCases;
            end
            for caseNum = caseNums
                disp(caseNum);
                cases = obj.paramCases(1).c;
                if caseNum ~= 1
                    cases = [cases,obj.paramCases(caseNum).c];
                end
                SIM21Analysis.plotGraphsByZ(cases);
                SIM21Analysis.plotEpsByZ(cases);
                SIM21Analysis.plotXeByZ(cases);
            end
        end
        
        
        function plotEpsGraphs(obj)
            for caseNum = obj.workCases
                disp(caseNum);
                SIM21Analysis.plotEpsByZ(obj.paramCases(caseNum).c);
            end
        end
        
        
        function plotXeGraphs(obj)
            for caseNum = obj.workCases
                disp(caseNum);
                SIM21Analysis.plotXeByZ(obj.paramCases(caseNum).c);
            end
        end


        function plotXeXHIGraphs(obj)
            for caseNum = obj.workCases
                disp(caseNum);
                SIM21Analysis.plotXeXHIByZ(obj.paramCases(caseNum).c);
            end
        end


        function plot4CasesZGraphs(obj,caseNums)
            outputPath = [obj.outputPath,'Graphs/Grouped4/'];
            if exist('caseNums','var')
                caseNums = intersect(caseNums,obj.workCases);
            else
                caseNums = obj.workCases;
            end
            for ind = 2:4:length(caseNums)
                disp([1,caseNums(ind:min(end,ind+3))]);
                cases = [obj.paramCases([1,caseNums(ind:min(end,ind+3))]).c];
                name = ['Cases',num2str(caseNums(ind)),'-',num2str(caseNums(min(end,ind+3)))];
                SIM21Analysis.plotGraphsByZ(cases,outputPath,name);
                SIM21Analysis.plotEpsByZ(cases,outputPath,name);
                SIM21Analysis.plotXeByZ(cases,outputPath,name);
            end
        end


        function plotSomething(obj)
            graphsXY = {{{'specialParams','minSlope','z'},{'specialParams','minSlope','T'}},...
                        {{'specialParams','maxSlope','z'},{'specialParams','maxSlope','T'}},...
                        {{'specialParams','minSlope','z'},{'specialParams','minSlope','slope'}},...
                        {{'specialParams','maxSlope','z'},{'specialParams','maxSlope','slope'}},...
                        {{'specialParams','minT21cm','z'},{'specialParams','minT21cm','T'}},...
                        {{'specialParams','maxT21cm','z'},{'specialParams','maxT21cm','T'}}};
            for xyFields = graphsXY
                obj.plotBasicSigScatter(xyFields{1}{1},xyFields{1}{2},[obj.outputPath,'Graphs/zT/']);
            end
        end

        function plotSomething4(obj)
            outputPath = [obj.outputPath,'Graphs/'];

            xfield = {'specialParams','minSlope','z'};
            yfield = {'specialParams','maxSlope','z'};
            figSettings = SIM21Analysis.initFigSettings('Max Slope - Min Slope - zMax(zMin)','',strjoin(xfield(2:end),'.'),strjoin(yfield(2:end),'.'));
            figSettings.plots = [figSettings.plots, obj.getLines(xfield,yfield,{'paramCases','c','fstar'})];
            obj.plotSigScatter(xfield,yfield,figSettings,[outputPath,'minMaxSlope.png']);

            xfield = {'paramCases','c','vc'};
            yfield = {'specialParams','maxSlope','z'};
            figSettings = SIM21Analysis.initFigSettings('Max Slope - vc,fstarVals','',strjoin(xfield(2:end),'.'),strjoin(yfield(2:end),'.'));
            figSettings.plots = [figSettings.plots, obj.getLines(xfield,yfield,{'paramCases','c','fstar'})];
            obj.plotSigScatter(xfield,yfield,figSettings,[outputPath,'MaxSlope - vc,fstar.png']);

            xfield = {'paramCases','c','vc'};
            yfield = {'tempParams','xHI75_25','zDiff'};
            figSettings = SIM21Analysis.initFigSettings('xHI75-xHI25 - vc','',strjoin(xfield(2:end),'.'),strjoin(yfield(2:end),'.'));
            obj.plotSigScatter(xfield,yfield,figSettings,[outputPath,'xHI75-xHI25 - vc.png']);

            xfield = {'paramCases','c','fstar'};
            yfield = {'specialParams','maxSlope','z'};
            figSettings = SIM21Analysis.initFigSettings('maxSlope - fstar','',strjoin(xfield(2:end),'.'),strjoin(yfield(2:end),'.'));
            obj.plotSigScatter(xfield,yfield,figSettings,[outputPath,'maxSlope-fstar.png']);
        end


        function findCorr(obj,corrLevel)
            corrParams = {{'specialParams','minT21cm','z'},{'specialParams','minT21cm','T'},...
                          {'specialParams','maxT21cm','z'},{'specialParams','maxT21cm','T'},...
                          {'specialParams','minSlope','z'},{'specialParams','minSlope','T'},{'specialParams','minSlope','slope'},...
                          {'specialParams','maxSlope','z'},{'specialParams','maxSlope','T'},{'specialParams','maxSlope','slope'},...
                          {'specialParams','minTK','z'},{'specialParams','minTK','T'},...
                          {'specialParams','xCross','z1'},...
                          {'specialParams','xHI75','z'},{'specialParams','xHI50','z'},{'specialParams','xHI25','z'},{'specialParams','xHI0','z'},...
                          {'specialParams','THT','z'},{'specialParams','THT','T'},...
                          {'paramCases','c','fstar'},{'paramCases','c','vbc'},{'paramCases','c','vc'},...
                          {'paramCases','c','fx'},{'paramCases','c','sed'},{'paramCases','atau'},...
                          {'tempParams','minMaxSlope','zDiff'},{'tempParams','xHI75_25','zDiff'}};
            for i = 1:length(corrParams)
                corrVectors(:,i) = obj.getField(obj.workCases,corrParams{i})';
            end

            corrs = triu(corr(corrVectors),1);
            indCorr = find(abs(corrs)>corrLevel);
            [row,col] = ind2sub(size(corrs),indCorr);

            outputPath = [obj.outputPath,'Graphs/Correlations/findcorr',num2str(corrLevel),'_'];
            for i = 1:length(row)
                xfield = corrParams{row(i)};
                yfield = corrParams{col(i)};
                xLabel = strjoin(xfield(2:end),'.');
                yLabel = strjoin(yfield(2:end),'.');
                figSettings = SIM21Analysis.initFigSettings([yLabel,' (',xLabel,') - ',num2str(corrs(indCorr(i)))],'',xLabel,yLabel);
                obj.plotSigScatter(xfield,yfield,figSettings,[outputPath,yLabel,'-',xLabel,'.png']);
            end
        end

        function plotBasicSigScatter(obj,xfield,yfield,outputPath,caseNums)
            [figSettings,outputName] = initSSFigSettings(obj,xfield,yfield,outputPath);
            if ~exist('caseNums','var')
                caseNums = [obj.regularCase,obj.smallVarCases,obj.largeVarCases];
            end
            plotSigScatter(obj,xfield,yfield,figSettings,outputName,caseNums);
        end


        function [figSettings,outputName] = initSSFigSettings(obj,xfield,yfield,outputPath,outputNamePrefix)
            xlabel = strjoin(xfield(2:end),'.');
            ylabel = strjoin(yfield(2:end),'.');
            figSettings = SIM21Analysis.initFigSettings([ylabel,'(',xlabel,')'],'',xlabel,ylabel);
            if ~exist('outputNamePrefix','var')
                outputNamePrefix = '';
            end
            outputName = [outputPath,outputNamePrefix,ylabel,'_',xlabel,'.png'];
        end


        function plotSigScatter(obj,xfield,yfield,figSettings,outputName,caseNums)
            caseTypes = {obj.regularCase,obj.smallVarCases,obj.largeVarCases,obj.otherCases};
            caseNames = {'regular','small variation','large variation','other'};
            %caseColors = {'r','g','b','k'};
            scatters = {};

            for caseTypeNum = 1:length(caseTypes)
                if exist('caseNums','var')
                    typeCaseNums = intersect(caseNums,caseTypes{caseTypeNum});
                else
                    typeCaseNums = caseTypes{caseTypeNum};
                end
                
                if ~ isempty(typeCaseNums)
                    x = obj.getField(typeCaseNums,xfield);
                    y = obj.getField(typeCaseNums,yfield);
                    scatters = [scatters,SIM21Analysis.plotScatter(cat(1,x,y),caseNames{caseTypeNum})];
                end
            end

            figSettings.plots = [scatters, figSettings.plots];
            SIM21Analysis.plotData(outputName,figSettings);
        end


        function lines = getLines(obj,xfield,yfield,lineField,caseNums)
            if ~ exist('caseNums','var')
                caseNums = obj.workCases;
            else
                caseNums = intersect(caseNums,obj.workCases);
            end

            lines = {};
            lineVals = obj.getField(caseNums,lineField);
            for lineVal = unique(lineVals)
                cases = caseNums(lineVals==lineVal);
                lines{end+1} = SIM21Analysis.plotLine(cat(1,obj.getField(cases,xfield),obj.getField(cases,yfield)),[lineField{end},' = ',num2str(lineVal)]);
            end
        end


        function vals = getField(obj,casesNum,fieldList)
            casesNum = obj.areRun(casesNum);
            vals = zeros(1,length(casesNum));
            for i = 1:length(vals)
                lst = getfield(obj,fieldList{1});
                val = lst(casesNum(i));

                for fieldID = 2:length(fieldList)
                    val = getfield(val,fieldList{fieldID});
                end

                if isempty(val)
                    vals(i) = NaN;
                else    
                    vals(i) = val;
                end
            end
        end
        

        function caseOutputPath = getCaseOutputPath(obj,caseNum)
            caseOutputPath = [obj.outputPath,'Case_',num2str(caseNum),'/'];
            [status,message,messageid] = mkdir(caseOutputPath);
        end
    end
end