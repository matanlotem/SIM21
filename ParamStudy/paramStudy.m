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
            obj.otherCases = [];
            
            obj.paramCases = obj.getCases(paramDataPath,cubeNum);
            obj.workCases = obj.areRun([obj.regularCase,obj.smallVarCases,obj.largeVarCases,obj.otherCases]);
            obj.specialParams = obj.getSpecialParams();
            obj.tempParams = obj.getTempParams();
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
            headers = ['General\t\t'   ,'minT21cm\t\t' ,'maxT21cm\t\t' ,'minSlope\t\t\t'      ,'maxSlope\t\t\t'      ,'minTK\t\t' ,'minTS\t\t\t'       ,'0 Crossing\t\t' ,'xHI percentage\t\t\t\t\t\t\t'                             ,'THT\t\t'   ,'final xHI\t\t'        ,'tau\t\t'              ,'\n',...
                       'Case\t','ID\t' ,'z\t','T\t'    ,'z\t','T\t'    ,'z\t','T\t','slope\t' ,'z\t','T\t','slope\t' ,'z\t','T\t' ,'z\t','T\t','xHI\t' ,'z\t','NOC\t'    ,'z80\t','z75\t','z50\t','z25\t','z0.05\t','z0.01\t','z0\t' ,'z\t','T\t' ,'planned\t','actual\t' ,'planned\t','actual\t' ,'\n'];
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
                          formatCell(obj.specialParams(caseNum).minTS.z)          ,'\t',...
                          formatCell(obj.specialParams(caseNum).minTS.T)          ,'\t',...
                          formatCell(obj.specialParams(caseNum).minTS.xHI)        ,'\t',...
                          formatCell(obj.specialParams(caseNum).xCross.z)         ,'\t',...
                          formatCell(length(obj.specialParams(caseNum).xCross.z)) ,'\t',...
                          formatCell(obj.specialParams(caseNum).xHI80.z)          ,'\t',...
                          formatCell(obj.specialParams(caseNum).xHI75.z)          ,'\t',...
                          formatCell(obj.specialParams(caseNum).xHI50.z)          ,'\t',...
                          formatCell(obj.specialParams(caseNum).xHI25.z)          ,'\t',...
                          formatCell(obj.specialParams(caseNum).xHI005.z)         ,'\t',...
                          formatCell(obj.specialParams(caseNum).xHI001.z)         ,'\t',...
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

        function plotAllT21cmByZGraphs(obj)
            outputPath = [obj.outputPath,'Graphs/zT/'];
            
            SIM21Analysis.plotT21cmByZ([obj.paramCases(obj.workCases).c],outputPath,'all');
            SIM21Analysis.plotT21cmByZ([obj.paramCases(intersect(obj.workCases,obj.smallVarCases)).c],outputPath,'smallVar');
            SIM21Analysis.plotT21cmByZ([obj.paramCases(intersect(obj.workCases,obj.largeVarCases)).c],outputPath,'largeVar');
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

        function plotSomething2(obj)
            graphsXY = {{{'specialParams','minSlope','z'},{'specialParams','minSlope','T'}},...
                        {{'specialParams','maxSlope','z'},{'specialParams','maxSlope','T'}},...
                        {{'specialParams','minT21cm','z'},{'specialParams','minT21cm','T'}},...
                        {{'specialParams','maxT21cm','z'},{'specialParams','maxT21cm','T'}}};
            
            prefix = 'd';
            outputPath = [obj.outputPath,'Graphs/zT/'];
            caseNums = [obj.regularCase,obj.smallVarCases,obj.largeVarCases];
            for ind = [1:4]
                xfield = graphsXY{ind}{1};
                yfield = graphsXY{ind}{2};
                [figSettings,outputName] = obj.initSSFigSettings(xfield,yfield,outputPath,[prefix,'_']);
                figSettings.plots = [figSettings.plots, obj.getLineByFieldVals(xfield,yfield,{{'paramCases','c','vc'},{'paramCases','c','fstar'}},[4.2,0.158],caseNums)];
                %figSettings.plots = [figSettings.plots, SIM21Analysis.plotLine([obj.getField([6,7,8,9,2,3,4,5],xfield); obj.getField([6,7,8,9,2,3,4,5],yfield)],'vc=4.2,fstar=0.158')];
                %figSettings.plots = [figSettings.plots, obj.getLineByFieldVals(xfield,yfield,{{'paramCases','c','vc'},{'paramCases','c','fstar'}},[35.5,0.015],caseNums)];
                figSettings.plots = [figSettings.plots, obj.getLineByFieldVals(xfield,yfield,{{'paramCases','c','vc'},{'paramCases','c','fstar'},{'paramCases','c','tau'}},[35.5,0.015,0.066],caseNums)];
                figSettings.plots = [figSettings.plots, obj.getLineByFieldVals(xfield,yfield,{{'paramCases','c','vc'},{'paramCases','c','fstar'},{'paramCases','c','tau'}},[35.5,0.015,0.082],caseNums)];
                %figSettings.plots = [figSettings.plots, SIM21Analysis.plotLine([obj.getField([30,32,26,28,29,27,33,31],xfield); obj.getField([30,32,26,28,29,27,33,31],yfield)],'vc=35.5,fstar=0.015')];
                figSettings.plots = [figSettings.plots, obj.getLineByFieldVals(xfield,yfield,{{'paramCases','c','fx'},{'paramCases','c','sed'}},[1.58,5],caseNums)];
                %figSettings.plots = [figSettings.plots, obj.getLineByFieldVals(xfield,yfield,{{'paramCases','c','fx'},{'paramCases','c','sed'}},[0.16,4],caseNums)];
                figSettings.plots = [figSettings.plots, obj.getLineByFieldVals(xfield,yfield,{{'paramCases','c','fx'},{'paramCases','c','sed'},{'paramCases','c','tau'}},[0.16,4,0.066],caseNums)]; 
                figSettings.plots = [figSettings.plots, obj.getLineByFieldVals(xfield,yfield,{{'paramCases','c','fx'},{'paramCases','c','sed'},{'paramCases','c','tau'}},[0.16,4,0.082  ],caseNums)];

                figSettings.plots = [figSettings.plots, obj.getLineByFieldVals(xfield,yfield,{{'paramCases','c','vc'},{'paramCases','c','fstar'}},[4.2,0.5],caseNums)]; 
                %figSettings.plots = [figSettings.plots, SIM21Analysis.plotLine([obj.getField([40,41,36,37,38,39],xfield); obj.getField([40,41,36,37,38,39],yfield)],'vc=4.2,fstar=0.5')];
                figSettings.plots = [figSettings.plots, obj.getLineByFieldVals(xfield,yfield,{{'paramCases','c','vc'},{'paramCases','c','fstar'}},[76.5,0.005],caseNums)];
                figSettings.plots = [figSettings.plots, obj.getLineByFieldVals(xfield,yfield,{{'paramCases','c','fx'},{'paramCases','c','sed'}},[10,2],caseNums)];
                %figSettings.plots = [figSettings.plots, obj.getLineByFieldVals(xfield,yfield,{{'paramCases','c','fx'},{'paramCases','c','sed'}},[0.1,6],caseNums)];
                figSettings.plots = [figSettings.plots, obj.getLineByFieldVals(xfield,yfield,{{'paramCases','c','fx'},{'paramCases','c','sed'},{'paramCases','c','tau'},{'paramCases','c','vc'}},[0.1,6,0.066,76.5],caseNums)];
                obj.plotSigScatter(xfield,yfield,figSettings,outputName,caseNums);
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


        function finalPlots(obj)
            %graphsXY = {{{'specialParams','minSlope','z'},{'specialParams','minSlope','T'}},...
            %            {{'specialParams','maxSlope','z'},{'specialParams','maxSlope','T'}},...
            %            {{'specialParams','minT21cm','z'},{'specialParams','minT21cm','T'}},...
            %            {{'specialParams','maxT21cm','z'},{'specialParams','maxT21cm','T'}}};
            
            outputPath = [obj.outputPath,'Graphs/Final/'];
            figSettings = SIM21Analysis.initFigSettings('','','$\nu\ [MHz]$','$T_b\ [mK]$');
            figSettings.legend = false;
            obj.plotFinalSigScatter({'specialParams','minT21cm','z'},{'specialParams','minT21cm','T'},figSettings,[outputPath,'Min21cm.png']);
            obj.plotFinalSigScatter({'specialParams','maxT21cm','z'},{'specialParams','maxT21cm','T'},figSettings,[outputPath,'Max21cm.png']);
            obj.plotFinalSigScatter({'specialParams','minSlope','z'},{'specialParams','minSlope','T'},figSettings,[outputPath,'MinSlopeT.png']);
            obj.plotFinalSigScatter({'specialParams','maxSlope','z'},{'specialParams','maxSlope','T'},figSettings,[outputPath,'MaxSlopeT.png']);
            figSettings = SIM21Analysis.initFigSettings('','','$\nu\ [MHz]$','$\frac{dT_b}{d\nu}\ [\frac{mK}{MHz}]$');
            figSettings.legend = false;
            obj.plotFinalSigScatter({'specialParams','minSlope','z'},{'specialParams','minSlope','slope'},figSettings,[outputPath,'MinSlopeS.png']);
            obj.plotFinalSigScatter({'specialParams','maxSlope','z'},{'specialParams','maxSlope','slope'},figSettings,[outputPath,'MaxSlopeS.png']);

            % Stack T21cm
            outputName = [outputPath,'T21cmStacked.png'];
            figSettings = SIM21Analysis.initFigSettings('','','$\nu\ [MHz]$','$T_b\ [mK]$');
            figSettings.legend = false;

            bigZetaCases = obj.areRun([58:65]);
            deltaZCases = obj.areRun([35:2:49]);
            tau93Cases = [34:2:48];
            badTSCases = obj.areRun(setdiff([40,41,48,49,56],[bigZetaCases,deltaZCases,tau93Cases]));
            smallVarCases = setdiff(obj.smallVarCases,[badTSCases,bigZetaCases,deltaZCases,tau93Cases]);
            largeVarCases = setdiff(obj.largeVarCases,[badTSCases,bigZetaCases,deltaZCases,tau93Cases]);

            T21cmDataReg = SIM21Analysis.getZData(obj.paramCases(1).c,'T21cm');
            T21cmDataReg = [SIM21Utils.T21cmWaveLength./(T21cmDataReg(1,:)+1); T21cmDataReg(2,:)];

            T21cmData56 = SIM21Analysis.getZData(obj.paramCases(56).c,'T21cm');
            figSettings.plots{end+1} = SIM21Analysis.plotLine([T21cmDataReg(1,:);T21cmData56(2,:)],'56',[0.9,0.9,0.9],'-',2);

            T21cmDatas = [];
            for caseNum = [deltaZCases,bigZetaCases];
                T21cmData = SIM21Analysis.getZData(obj.paramCases(caseNum).c,'T21cm');
                T21cmDatas = [T21cmDatas; T21cmData(2,:)];
            end
            figSettings.plots{end+1} = SIM21Analysis.plotFill([T21cmDataReg(1,:); min(T21cmDatas)],[T21cmDataReg(1,:); max(T21cmDatas)],'bigZeta and delta Z',[0.9,0.9,0.9]);

            T21cmDatas = [];
            for caseNum = largeVarCases;
                T21cmData = SIM21Analysis.getZData(obj.paramCases(caseNum).c,'T21cm');
                T21cmDatas = [T21cmDatas; T21cmData(2,:)];
            end
            figSettings.plots{end+1} = SIM21Analysis.plotFill([T21cmDataReg(1,:); min(T21cmDatas)],[T21cmDataReg(1,:); max(T21cmDatas)],'Large Variations',[0.7,0.7,0.7]);

            T21cmDatas = [];
            for caseNum = smallVarCases
                T21cmData = SIM21Analysis.getZData(obj.paramCases(caseNum).c,'T21cm');
                T21cmDatas = [T21cmDatas; T21cmData(2,:)];
            end
            figSettings.plots{end+1} = SIM21Analysis.plotFill([T21cmDataReg(1,:); min(T21cmDatas)],[T21cmDataReg(1,:); max(T21cmDatas)],'Small Variations',[0.5,0.5,0.5]);

            figSettings.plots{end+1} = SIM21Analysis.plotLine(T21cmDataReg,'Reg','k','-',2);

            f = figure();
            figSettings.xTick = [0:30:300];
            SIM21Analysis.plotPlot(f,figSettings);
            axnu = gca;
            axz = axes('Position',axnu.Position,...
                       'XAxisLocation','top',...
                       'Color','none');
            axz.YTick = [];
            axnuXLim = axnu.XLim;
            gca;
            xlim(axnuXLim);
            maxZ1 = SIM21Utils.T21cmWaveLength/axnuXLim(1);
            minZ1 = SIM21Utils.T21cmWaveLength/axnuXLim(2);
            xTickVals = flip(floor(ceil(minZ1).*(1.5.^[0:8])));
            axz.XTick = SIM21Utils.T21cmWaveLength./xTickVals;
            axz.XTickLabel = num2str(xTickVals');
            xlabel('$1+z$','FontSize',12,'Interpreter','LaTex');
            saveas(f,outputName);
            saveas(f,outputName(1:end-4),'epsc');

        end

        function plotFinalSigScatter(obj,xfield,yfield,figSettings,outputName)
            bigZetaCases = obj.areRun([58:65]);
            deltaZCases = obj.areRun([35:2:49]);
            tau93Cases = [34:2:48];
            badTSCases = obj.areRun(setdiff([40,41,48,49,56],[bigZetaCases,deltaZCases,tau93Cases]));
            smallVarCases = setdiff(obj.smallVarCases,[badTSCases,bigZetaCases,deltaZCases,tau93Cases]);
            largeVarCases = setdiff(obj.largeVarCases,[badTSCases,bigZetaCases,deltaZCases,tau93Cases]);

            caseTypes = {obj.regularCase,smallVarCases,largeVarCases,badTSCases,deltaZCases,bigZetaCases};
            caseNames = {'regular','small variations','large variations','bad Ts','\Delta z','bad zeta'};

            %caseColors = {[0.3,0.4,0.9],[0.7,0.15,0.15],[0.9,0.75,0],[0.79,0.9,0.79],[0.87,0.73,0.9],[0.8,0.8,0.8]};
            caseColors = {[0.3,0.4,0.9],[0.7,0.15,0.15],[0.9,0.75,0],[0.8,0.8,0.8],[0.8,0.8,0.8],[0.8,0.8,0.8]};
            caseShapes = {'o','o','o','p','s','d'};
            scatters = {};

            for caseTypeNum = 1:length(caseTypes)
                typeCaseNums = caseTypes{caseTypeNum};
                
                if ~ isempty(typeCaseNums)
                    x = obj.getField(typeCaseNums,xfield);
                    x = SIM21Utils.T21cmWaveLength./(1+x);
                    y = obj.getField(typeCaseNums,yfield);
                    scatters{end+1} = SIM21Analysis.plotScatter([x;y],caseNames{caseTypeNum},caseColors{caseTypeNum},caseShapes{caseTypeNum});
                end
            end

            figSettings.plots = [scatters, figSettings.plots];
            figSettings.xTick = [20:20:400];
            f = figure();
            SIM21Analysis.plotPlot(f,figSettings);
            axnu = gca;
            axz = axes('Position',axnu.Position,...
                       'XAxisLocation','top',...
                       'Color','none');
            axz.YTick = [];
            axnuXLim = axnu.XLim;
            gca;
            xlim(axnuXLim);
            maxZ1 = SIM21Utils.T21cmWaveLength/axnuXLim(1);
            minZ1 = SIM21Utils.T21cmWaveLength/axnuXLim(2);
            xTickVals = flip(floor(ceil(minZ1).*(1.3.^[0:5])));
            axz.XTick = SIM21Utils.T21cmWaveLength./xTickVals;
            axz.XTickLabel = num2str(xTickVals');
            xlabel('$1+z$','FontSize',12,'Interpreter','LaTex');
            saveas(f,outputName);
            saveas(f,outputName(1:end-4),'epsc');
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

        function line = getLineByFieldVals(obj,xfield,yfield,lineFields,lineVals,caseNums)
            cases = obj.getCasesByFieldVals(lineFields,lineVals,caseNums);
            name = [];
            for ind = 1:length(lineFields)
                name=[name,',',lineFields{ind}{end},'=',num2str(lineVals(ind))];
            end
            name = name(2:end);
            if isempty(name)
                name = '';
            end

            line = SIM21Analysis.plotLine(sortrows([obj.getField(cases,xfield);obj.getField(cases,yfield)]')',name);
            %line = SIM21Analysis.plotLine(cat(1,obj.getField(cases,xfield),obj.getField(cases,yfield)),name);
        end

        function cases = getCasesByFieldVals(obj,fields,vals,caseNums)
            if ~ exist('caseNums','var')
                caseNums = obj.workCases;
            else
                caseNums = intersect(caseNums,obj.workCases);
            end

            name = [];
            cases = caseNums(:)';
            for ind = 1:length(fields)
                cases = intersect(cases,caseNums(obj.getField(caseNums,fields{ind})==vals(ind)));
            end
        end


        function vals = getField(obj,caseNums,fieldList)
            caseNums = obj.areRun(caseNums);
            vals = zeros(1,length(caseNums));
            for i = 1:length(vals)
                lst = getfield(obj,fieldList{1});
                val = lst(caseNums(i));

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