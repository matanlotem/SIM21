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
        function obj = paramStudy(ps)
            addpath('/a/home/cc/tree/taucc/students/physics/matanlotem/Work/SIM21/lib/');
            if ~exist('ps','var')
                ps = 2;
            end

            switch ps
                case 1
                    paramDataPath = 'paramStudy.xlsx';
                    cubeNum = 9;
                    obj.pathExt = 'PS_';
                    obj.outputPath = '/scratch300/matanlotem/ParamStudy/';
                    obj.regularCase = [1];
                    obj.smallVarCases = [2:33];
                    obj.largeVarCases = [37:2:49 50:55 57:61 63 65];
                    obj.otherCases = [];
                    
                    obj.paramCases = obj.getCases(paramDataPath,cubeNum);
                    obj.workCases = obj.areRun([obj.regularCase,obj.smallVarCases,obj.largeVarCases,obj.otherCases]);
                    obj.specialParams = obj.getSpecialParams();
                    obj.tempParams = obj.getTempParams();
                case 2
                    paramDataPath = 'paramStudy2.xlsx';
                    cubeNum = 9;
                    obj.pathExt = 'GRID_';
                    obj.outputPath = '/scratch300/matanlotem/GRID/';
                    % obj.regularCase = [1];
                    % obj.smallVarCases = [2:33];
                    % obj.largeVarCases = [37:2:49 50:65];
                    % obj.otherCases = [];
                    
                    obj.paramCases = obj.getCases(paramDataPath,cubeNum);
                    % obj.workCases = obj.areRun([obj.regularCase,obj.smallVarCases,obj.largeVarCases,obj.otherCases]);
                    % obj.specialParams = obj.getSpecialParams();
                    % obj.tempParams = obj.getTempParams();
            end
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
                tempParams(caseNum).vc3fstar = obj.paramCases(caseNum).c.vc^3 * obj.paramCases(caseNum).c.fstar;
                tempParams(caseNum).fmaxT21cmFstarFcoll = obj.paramCases(caseNum).c.fstar * obj.specialParams(caseNum).fmaxT21cm.fcoll;
                tempParams(caseNum).minT21cmFstarFcollFx = obj.paramCases(caseNum).c.fstar * obj.specialParams(caseNum).minT21cm.fcoll * obj.paramCases(caseNum).c.fx;
                tempParams(caseNum).minT21cmFstarFcoll = obj.paramCases(caseNum).c.fstar * obj.specialParams(caseNum).minT21cm.fcoll;
                tempParams(caseNum).maxT21cmZetaFcoll = obj.paramCases(caseNum).c.zeta * obj.specialParams(caseNum).maxT21cm.fcoll;
                tempParams(caseNum).maxT21cmZetaFcoll_div_tau = obj.paramCases(caseNum).c.zeta * obj.specialParams(caseNum).maxT21cm.fcoll / obj.paramCases(caseNum).atau;
                tempParams(caseNum).maxT21cmZetaFcoll_div_tau2 = obj.paramCases(caseNum).c.zeta * obj.specialParams(caseNum).maxT21cm.fcoll / obj.paramCases(caseNum).atau^2;
                tempParams(caseNum).maxT21cmZetaFcoll_div_tau3 = obj.paramCases(caseNum).c.zeta * obj.specialParams(caseNum).maxT21cm.fcoll / obj.paramCases(caseNum).atau^3;
                tempParams(caseNum).maxT21cmZetaFcoll_div_tau4 = obj.paramCases(caseNum).c.zeta * obj.specialParams(caseNum).maxT21cm.fcoll / obj.paramCases(caseNum).atau^4;
                tempParams(caseNum).maxT21cmZetaFcoll_div_tau10 = obj.paramCases(caseNum).c.zeta * obj.specialParams(caseNum).maxT21cm.fcoll / obj.paramCases(caseNum).atau^10;
                tempParams(caseNum).div_tau = 1 / obj.paramCases(caseNum).atau;
            end
        end
        

        function specialParamsTable(obj)
            headers = ['General\t\t'   ,'minTb\t\t\t'         ,'maxTb\t\t\t'         ,'firstMaxTb\t\t\t'    ,'minTbSlope\t\t\t'    ,'maxTbSlope\t\t\t'    ,'minTK\t\t' ,'minTS\t\t\t'       ,'0 Crossing\t\t' ,'xHI percentage\t\t\t\t\t\t\t\t'                                   ,'THT\t\t'   ,'Tb=0\t' ,'final xHI\t\t'        ,'tau\t\t'              ,'\n',...
                       'Case\t','ID\t' ,'z\t','T\t','fcoll\t' ,'z\t','T\t','fcoll\t' ,'z\t','T\t','fcoll\t' ,'z\t','T\t','slope\t' ,'z\t','T\t','slope\t' ,'z\t','T\t' ,'z\t','T\t','xHI\t' ,'z\t','NOC\t'    ,'z95\t','z80\t','z75\t','z50\t','z25\t','z0.05\t','z0.01\t','z0\t' ,'z\t','T\t' ,'z\t'    ,'planned\t','actual\t' ,'planned\t','actual\t' ,'\n'];
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
                          formatCell(obj.specialParams(caseNum).minT21cm.fcoll)   ,'\t',...
                          formatCell(obj.specialParams(caseNum).maxT21cm.z)       ,'\t',...
                          formatCell(obj.specialParams(caseNum).maxT21cm.T)       ,'\t',...
                          formatCell(obj.specialParams(caseNum).maxT21cm.fcoll)   ,'\t',...
                          formatCell(obj.specialParams(caseNum).fmaxT21cm.z)      ,'\t',...
                          formatCell(obj.specialParams(caseNum).fmaxT21cm.T)      ,'\t',...
                          formatCell(obj.specialParams(caseNum).fmaxT21cm.fcoll)  ,'\t',...
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
                          formatCell(obj.specialParams(caseNum).xHI95.z)          ,'\t',...
                          formatCell(obj.specialParams(caseNum).xHI80.z)          ,'\t',...
                          formatCell(obj.specialParams(caseNum).xHI75.z)          ,'\t',...
                          formatCell(obj.specialParams(caseNum).xHI50.z)          ,'\t',...
                          formatCell(obj.specialParams(caseNum).xHI25.z)          ,'\t',...
                          formatCell(obj.specialParams(caseNum).xHI005.z)         ,'\t',...
                          formatCell(obj.specialParams(caseNum).xHI001.z)         ,'\t',...
                          formatCell(obj.specialParams(caseNum).xHI0.z)           ,'\t',...
                          formatCell(obj.specialParams(caseNum).THT.z)            ,'\t',...
                          formatCell(obj.specialParams(caseNum).THT.T)            ,'\t',...
                          formatCell(obj.specialParams(caseNum).T21cm0.z)         ,'\t',...
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
            %%bigZetaCases = obj.areRun([58:65]);
            %%deltaZCases = obj.areRun([35:2:49]);
            tau93Cases = [34:2:48];
            %%badTSCases = obj.areRun(setdiff([40,41,48,49,56],[bigZetaCases,deltaZCases,tau93Cases]));
            %%smallVarCases = setdiff(obj.smallVarCases,[badTSCases,bigZetaCases,deltaZCases,tau93Cases]);
            %%largeVarCases = setdiff(obj.largeVarCases,[badTSCases,bigZetaCases,deltaZCases,tau93Cases]);

            %%caseTypes = {obj.regularCase,smallVarCases,largeVarCases,badTSCases,deltaZCases,bigZetaCases};
            %%caseNames = {'regular','small variations','large variations','bad Ts','\Delta z','bad zeta'};
            %%caseColors = {[0.3,0.4,0.9],[0.7,0.15,0.15],[0.9,0.75,0],[0.8,0.8,0.8],[0.8,0.8,0.8],[0.8,0.8,0.8]};
            %%caseShapes = {'o','o','o','p','s','d'};

            caseTypes = {obj.regularCase,obj.smallVarCases,obj.largeVarCases};
            caseNames = {'regular','small variations','large variations'};
            caseColors = {[0.3,0.4,0.9],[0.7,0.15,0.15],[0.9,0.75,0]};
            caseShapes = {'o','o','o'};

            tau111Case = [35];
            tau66Cases = setdiff(obj.getCasesByFieldVals({{'paramCases','c','tau'}},[0.066],obj.workCases),[tau111Case,tau93Cases]);
            tau82Cases = setdiff(obj.getCasesByFieldVals({{'paramCases','c','tau'}},[0.082],obj.workCases),[tau111Case,tau93Cases]);
            tau98Cases = setdiff(obj.getCasesByFieldVals({{'paramCases','c','tau'}},[0.098],obj.workCases),[tau111Case,tau93Cases]);

            caseTypes = {obj.regularCase,intersect(obj.smallVarCases,tau66Cases),intersect(obj.smallVarCases,tau82Cases),intersect(obj.largeVarCases,tau66Cases),intersect(obj.largeVarCases,tau98Cases)};
            caseTypesTmp = {obj.regularCase,intersect(obj.smallVarCases,tau66Cases),intersect(obj.smallVarCases,tau82Cases),intersect(obj.largeVarCases,tau66Cases),setdiff(intersect(obj.largeVarCases,tau98Cases),[37,41,45])};
            caseNames = {'regular','small variations - tau 66','small variations - tau 82','large variations - tau 66','large variations - tau 98'};
            caseColors = {[0.3,0.4,0.9],[0.7,0.15,0.15],[0.7,0.15,0.15],[0.9,0.75,0],[0.9,0.75,0]};
            %caseShapes = {'o','o','s','o','d'};
            caseShapes = {'o','o','o','o','o'};
            caseFilled = {0,0,1,0,1};

            tauCaseTypes = {tau66Cases,tau82Cases,tau98Cases};
            tauCaseNames = {'tau=0.066','tau=0.082','tau=0.098'};
            tauCaseColors = {[0.3,0.4,0.9],[0.7,0.15,0.15],[0.9,0.75,0]};%,[0.2,0.8,0.6]};
            tauCaseShapes = {'o','o','o'};
            tauCaseFilled = {0,0,0};
            
            outputPath = [obj.outputPath,'Graphs/Final/'];

            % Points of interest by redshift - Tb
            figSettings = SIM21Analysis.initFigSettings('','','$\nu\ [MHz]$','$T_b\ [mK]$',0);
            obj.plotFinalSigScatter({'specialParams','minT21cm','z'},{'specialParams','minT21cm','T'},figSettings,[outputPath,'MinTb.png'],caseTypes,caseNames,caseColors,caseShapes,caseFilled,1);
            %%obj.plotFinalSigScatter({'specialParams','minT21cm','z'},{'specialParams','minT21cm','T'},figSettings,[outputPath,'Tau_MinT21cm.png'],tauCaseTypes,tauCaseNames,tauCaseColors,tauCaseShapes,tauCaseFilled,1);
            obj.plotFinalSigScatter({'specialParams','maxT21cm','z'},{'specialParams','maxT21cm','T'},figSettings,[outputPath,'MaxTb.png'],caseTypes,caseNames,caseColors,caseShapes,caseFilled,1);
            %%obj.plotFinalSigScatter({'specialParams','maxT21cm','z'},{'specialParams','maxT21cm','T'},figSettings,[outputPath,'Tau_Max21cm.png'],tauCaseTypes,tauCaseNames,tauCaseColors,tauCaseShapes,tauCaseFilled,1);
            %$obj.plotFinalSigScatter({'specialParams','maxSlope','z'},{'specialParams','maxSlope','T'},figSettings,[outputPath,'MaxSlopeT.png'],caseTypes,caseNames,caseColors,caseShapes,1);
            %$obj.plotFinalSigScatter({'specialParams','maxSlope','z'},{'specialParams','maxSlope','T'},figSettings,[outputPath,'Tau_MaxSlopeT.png'],tauCaseTypes,tauCaseNames,tauCaseColors,tauCaseShapes,1);
            obj.plotFinalSigScatter({'specialParams','fmaxT21cm','z'},{'specialParams','fmaxT21cm','T'},figSettings,[outputPath,'FirstMaxTb.png'],caseTypes,caseNames,caseColors,caseShapes,caseFilled,1);
            %%obj.plotFinalSigScatter({'specialParams','fmaxT21cm','z'},{'specialParams','fmaxT21cm','T'},figSettings,[outputPath,'Tau_FirstMax21cm.png'],tauCaseTypes,tauCaseNames,tauCaseColors,tauCaseShapes,tauCaseFilled,1);
            %%% add vc-fstar lines
            %%c = [obj.paramCases(obj.workCases).c];
            %%for vc = unique([c.vc])
            %%    for fstar = unique([c.fstar])
            %%        l = obj.getLineByFieldVals({'specialParams','minSlope','z'},{'specialParams','minSlope','T'},{{'paramCases','c','vc'},{'paramCases','c','fstar'}},[vc,fstar],obj.workCases);
            %%        if ~ isempty(l.x)
            %%            l.x = l.x+1;
            %%            figSettings.plots = [figSettings.plots, l];
            %%        end
            %%    end
            %%end
            %%obj.plotFinalSigScatter({'specialParams','minSlope','z'},{'specialParams','minSlope','T'},figSettings,[outputPath,'MinSlopeT.png'],caseTypes,caseNames,caseColors,caseShapes,1);
            %%obj.plotFinalSigScatter({'specialParams','minSlope','z'},{'specialParams','minSlope','T'},figSettings,[outputPath,'Tau_MinSlopeT.png'],tauCaseTypes,tauCaseNames,tauCaseColors,tauCaseShapes,1);

            %%% Points of interest by redshift - slope
            %%figSettings = SIM21Analysis.initFigSettings('','','$\nu\ [MHz]$','$\frac{dT_b}{d\nu}\ [\frac{mK}{MHz}]$');
            %%figSettings.legend = false;
            %%obj.plotFinalSigScatter({'specialParams','maxSlope','z'},{'specialParams','maxSlope','slope'},figSettings,[outputPath,'MaxSlopeS.png'],caseTypes,caseNames,caseColors,caseShapes,1);
            %%obj.plotFinalSigScatter({'specialParams','maxSlope','z'},{'specialParams','maxSlope','slope'},figSettings,[outputPath,'Tau_MaxSlopeS.png'],tauCaseTypes,tauCaseNames,tauCaseColors,tauCaseShapes,1);
            %%% add vc-fstar lines
            %%c = [obj.paramCases(obj.workCases).c];
            %%for vc = unique([c.vc])
            %%    for fstar = unique([c.fstar])
            %%        l = obj.getLineByFieldVals({'specialParams','minSlope','z'},{'specialParams','minSlope','slope'},{{'paramCases','c','vc'},{'paramCases','c','fstar'}},[vc,fstar],obj.workCases);
            %%        if ~ isempty(l.x)
            %%            l.x = l.x+1;
            %%            figSettings.plots = [figSettings.plots, l];
            %%        end
            %%    end
            %%end
            %%obj.plotFinalSigScatter({'specialParams','minSlope','z'},{'specialParams','minSlope','slope'},figSettings,[outputPath,'MinSlopeS.png'],caseTypes,caseNames,caseColors,caseShapes,1);
            %%obj.plotFinalSigScatter({'specialParams','minSlope','z'},{'specialParams','minSlope','slope'},figSettings,[outputPath,'Tau_MinSlopeS.png'],tauCaseTypes,tauCaseNames,tauCaseColors,tauCaseShapes,1);
            
            % fcoll at Tb extremum
            figSettings = SIM21Analysis.initFigSettings('','','$\nu\ [MHz]$','$fstar * fcoll * fx$',0);
            figSettings.log = 'y';
            figSettings.yTick = 10.^[-10:0];
            obj.plotFinalSigScatter({'specialParams','minT21cm','z'},{'tempParams','minT21cmFstarFcollFx'},figSettings,[outputPath,'minTbZ-FstarFcollFx.png'],caseTypes,caseNames,caseColors,caseShapes,caseFilled,1);
            figSettings = SIM21Analysis.initFigSettings('','','$T_b\ [mK]$','$fstar * fcoll * fx$',0);
            figSettings.log = 'y';
            figSettings.xTick = [-300:50:100];
            figSettings.yTick = 10.^[-10:0];
            obj.plotFinalSigScatter({'specialParams','minT21cm','T'},{'tempParams','minT21cmFstarFcollFx'},figSettings,[outputPath,'minTbT-FstarFcollFx.png'],caseTypes,caseNames,caseColors,caseShapes,caseFilled);

            figSettings = SIM21Analysis.initFigSettings('','','$\nu\ [MHz]$','$fstar * fcoll$',0);
            figSettings.log = 'y';
            figSettings.yTick = 10.^[-10:0];
            obj.plotFinalSigScatter({'specialParams','minT21cm','z'},{'tempParams','minT21cmFstarFcoll'},figSettings,[outputPath,'minTbZ-FstarFcoll.png'],caseTypes,caseNames,caseColors,caseShapes,caseFilled,1);
            figSettings = SIM21Analysis.initFigSettings('','','$T_b\ [mK]$','$fstar * fcoll$',0);
            figSettings.log = 'y';
            figSettings.xTick = [-300:50:100];
            figSettings.yTick = 10.^[-10:0];
            obj.plotFinalSigScatter({'specialParams','minT21cm','T'},{'tempParams','minT21cmFstarFcoll'},figSettings,[outputPath,'minTbT-FstarFcoll.png'],caseTypes,caseNames,caseColors,caseShapes,caseFilled);


            figSettings = SIM21Analysis.initFigSettings('','','$\nu\ [MHz]$','$fstar * fcoll$',0);
            figSettings.log = 'y';
            figSettings.yTick = 10.^[-15:0];
            obj.plotFinalSigScatter({'specialParams','fmaxT21cm','z'},{'tempParams','fmaxT21cmFstarFcoll'},figSettings,[outputPath,'fmaxTbZ-FstarFcoll.png'],caseTypes,caseNames,caseColors,caseShapes,caseFilled,1);
            figSettings = SIM21Analysis.initFigSettings('','','$T_b\ [mK]$','$fstar * fcoll$',0);
            figSettings.xTick = [-40:5:40];
            figSettings.log = 'y';
            figSettings.yTick = 10.^[-15:0];
            obj.plotFinalSigScatter({'specialParams','fmaxT21cm','T'},{'tempParams','fmaxT21cmFstarFcoll'},figSettings,[outputPath,'fmaxTbT-FstarFcoll.png'],caseTypes,caseNames,caseColors,caseShapes,caseFilled);

            figSettings = SIM21Analysis.initFigSettings('','','$\nu\ [MHz]$','$zeta * fcoll$',0);
            obj.plotFinalSigScatter({'specialParams','maxT21cm','z'},{'tempParams','maxT21cmZetaFcoll'},figSettings,[outputPath,'maxTbZ-ZetaFcoll.png'],caseTypes,caseNames,caseColors,caseShapes,caseFilled,1);
            figSettings = SIM21Analysis.initFigSettings('','','$T_b\ [mK]$','$zeta * fcoll$',0);
            figSettings.xTick = [-40:5:40];
            obj.plotFinalSigScatter({'specialParams','maxT21cm','T'},{'tempParams','maxT21cmZetaFcoll'},figSettings,[outputPath,'maxTbT-ZetaFcoll.png'],caseTypes,caseNames,caseColors,caseShapes,caseFilled);

            % / tau
            figSettings = SIM21Analysis.initFigSettings('','','$\nu\ [MHz]$','$zeta * fcoll / \tau$',0);
            obj.plotFinalSigScatter({'specialParams','maxT21cm','z'},{'tempParams','maxT21cmZetaFcoll_div_tau'},figSettings,[outputPath,'maxTbZ-ZetaFcoll_div_tau.png'],caseTypes,caseNames,caseColors,caseShapes,caseFilled,1);
            figSettings = SIM21Analysis.initFigSettings('','','$T_b\ [mK]$','$zeta * fcoll / \tau$',0);
            figSettings.xTick = [-40:5:40];
            obj.plotFinalSigScatter({'specialParams','maxT21cm','T'},{'tempParams','maxT21cmZetaFcoll_div_tau'},figSettings,[outputPath,'maxTbT-ZetaFcoll_div_tau.png'],caseTypes,caseNames,caseColors,caseShapes,caseFilled);


            % / tau^2
            figSettings = SIM21Analysis.initFigSettings('','','$\nu\ [MHz]$','$zeta * fcoll / \tau^2$',0);
            obj.plotFinalSigScatter({'specialParams','maxT21cm','z'},{'tempParams','maxT21cmZetaFcoll_div_tau2'},figSettings,[outputPath,'maxTbZ-ZetaFcoll_div_tau2.png'],caseTypes,caseNames,caseColors,caseShapes,caseFilled,1);
            figSettings = SIM21Analysis.initFigSettings('','','$T_b\ [mK]$','$zeta * fcoll / \tau^2$',0);
            figSettings.xTick = [-40:5:40];
            obj.plotFinalSigScatter({'specialParams','maxT21cm','T'},{'tempParams','maxT21cmZetaFcoll_div_tau2'},figSettings,[outputPath,'maxTbT-ZetaFcoll_div_tau2.png'],caseTypes,caseNames,caseColors,caseShapes,caseFilled);

            % / tau^3
            figSettings = SIM21Analysis.initFigSettings('','','$\nu\ [MHz]$','$zeta * fcoll / \tau^3$',0);
            obj.plotFinalSigScatter({'specialParams','maxT21cm','z'},{'tempParams','maxT21cmZetaFcoll_div_tau3'},figSettings,[outputPath,'maxTbZ-ZetaFcoll_div_tau3.png'],caseTypes,caseNames,caseColors,caseShapes,caseFilled,1);
            figSettings = SIM21Analysis.initFigSettings('','','$T_b\ [mK]$','$zeta * fcoll / \tau^3$',0);
            figSettings.xTick = [-40:5:40];
            obj.plotFinalSigScatter({'specialParams','maxT21cm','T'},{'tempParams','maxT21cmZetaFcoll_div_tau3'},figSettings,[outputPath,'maxTbT-ZetaFcoll_div_tau3.png'],caseTypes,caseNames,caseColors,caseShapes,caseFilled);

            % / tau^4
            figSettings = SIM21Analysis.initFigSettings('','','$\nu\ [MHz]$','$zeta * fcoll / \tau^4$',0);
            obj.plotFinalSigScatter({'specialParams','maxT21cm','z'},{'tempParams','maxT21cmZetaFcoll_div_tau4'},figSettings,[outputPath,'maxTbZ-ZetaFcoll_div_tau4.png'],caseTypes,caseNames,caseColors,caseShapes,caseFilled,1);
            figSettings = SIM21Analysis.initFigSettings('','','$T_b\ [mK]$','$zeta * fcoll / \tau^4$',0);
            figSettings.xTick = [-40:5:40];
            obj.plotFinalSigScatter({'specialParams','maxT21cm','T'},{'tempParams','maxT21cmZetaFcoll_div_tau4'},figSettings,[outputPath,'maxTbT-ZetaFcoll_div_tau4.png'],caseTypes,caseNames,caseColors,caseShapes,caseFilled);


            % / tau^10
            figSettings = SIM21Analysis.initFigSettings('','','$\nu\ [MHz]$','$zeta * fcoll / \tau^10$',0);
            obj.plotFinalSigScatter({'specialParams','maxT21cm','z'},{'tempParams','maxT21cmZetaFcoll_div_tau10'},figSettings,[outputPath,'maxTbZ-ZetaFcoll_div_tau10.png'],caseTypes,caseNames,caseColors,caseShapes,caseFilled,1);
            figSettings = SIM21Analysis.initFigSettings('','','$T_b\ [mK]$','$zeta * fcoll / \tau^10$',0);
            figSettings.xTick = [-40:5:40];
            obj.plotFinalSigScatter({'specialParams','maxT21cm','T'},{'tempParams','maxT21cmZetaFcoll_div_tau10'},figSettings,[outputPath,'maxTbT-ZetaFcoll_div_tau10.png'],caseTypes,caseNames,caseColors,caseShapes,caseFilled);

            % 1 / tau
            figSettings = SIM21Analysis.initFigSettings('','','$\nu\ [MHz]$','$1 / \tau$',0);
            obj.plotFinalSigScatter({'specialParams','maxT21cm','z'},{'tempParams','div_tau'},figSettings,[outputPath,'maxTbZ-1_div_tau.png'],caseTypes,caseNames,caseColors,caseShapes,caseFilled,1);
            figSettings = SIM21Analysis.initFigSettings('','','$T_b\ [mK]$','$1 / \tau$',0);
            figSettings.xTick = [-40:5:40];
            obj.plotFinalSigScatter({'specialParams','maxT21cm','T'},{'tempParams','div_tau'},figSettings,[outputPath,'maxTbT-1_div_tau.png'],caseTypes,caseNames,caseColors,caseShapes,caseFilled);

            % Reionization points of interest by redshft to Tb=0
            %%figSettings = SIM21Analysis.initFigSettings('','','$T_b = 0\ \ (z+1)$','',0);
            %%figSettings.yLabel = '$xHI = 0.95\ \ (z+1)$';
            %%obj.plotFinalSigScatter({'specialParams','T21cm0','z'},{'specialParams','xHI95','z'},figSettings,[outputPath,'THTxHI95.png'],caseTypes,caseNames,caseColors,caseShapes,caseFilled);
            %%obj.plotFinalSigScatter({'specialParams','T21cm0','z'},{'specialParams','xHI95','z'},figSettings,[outputPath,'Tau_THTxHI95.png'],tauCaseTypes,tauCaseNames,tauCaseColors,tauCaseShapes,tauCaseFilled);
            %%figSettings.yLabel = '$xHI = 0.8\ \ (z+1)$';
            %%obj.plotFinalSigScatter({'specialParams','T21cm0','z'},{'specialParams','xHI80','z'},figSettings,[outputPath,'THTxHI80.png'],caseTypes,caseNames,caseColors,caseShapes,caseFilled);
            %%obj.plotFinalSigScatter({'specialParams','T21cm0','z'},{'specialParams','xHI80','z'},figSettings,[outputPath,'Tau_THTxHI80.png'],tauCaseTypes,tauCaseNames,tauCaseColors,tauCaseShapes,tauCaseFilled);
            %%figSettings.yLabel = '$xHI = 0.5\ \ (z+1)$';
            %%obj.plotFinalSigScatter({'specialParams','T21cm0','z'},{'specialParams','xHI50','z'},figSettings,[outputPath,'THTxHI50.png'],caseTypes,caseNames,caseColors,caseShapes,caseFilled);
            %%obj.plotFinalSigScatter({'specialParams','T21cm0','z'},{'specialParams','xHI50','z'},figSettings,[outputPath,'Tau_THTxHI50.png'],tauCaseTypes,tauCaseNames,tauCaseColors,tauCaseShapes,tauCaseFilled);

            %%figSettings = SIM21Analysis.initFigSettings('','','$1+z$','${vc}^3*fstar$');
            %%figSettings.log = 'xy';
            %%obj.plotFinalSigScatter({'specialParams','fmaxT21cm','z'},{'tempParams','vc3fstar'},figSettings,[outputPath,'vc3fstar_z.png'],caseTypes,caseNames,caseColors,caseShapes,caseFilled);
            

            % Stack T21cm
            %%stackCaseTypes = {[56],[deltaZCases,bigZetaCases],largeVarCases,smallVarCases,[1]};
            %%stackCaseNames = {'56','bigZeta and delta Z','Large Variations','Small Variations','Regular'};
            %%stackCaseColors = {[0.9,0.9,0.9],[0,1,0],[1,0,0],[0,0,1],[0,0,0]};
            stackCaseTypes = {[56],setdiff(obj.largeVarCases,[56]),obj.smallVarCases,obj.regularCase};
            stackCaseNames = {'56','Large Variations','Small Variations','Regular'};
            stackCaseColors = {[0.9,0.9,0.9],[1,0,0],[0,0,1],[0,0,0]};

            figSettings = SIM21Analysis.initFigSettings('','','$\nu\ [MHz]$','$T_b\ [mK]$',0);
            obj.plotStackedFills(figSettings,[outputPath,'TbStacked.png'],stackCaseTypes,stackCaseNames,stackCaseColors);
            obj.plotStackedFills(figSettings,[outputPath,'Tau_TbStacked.png'],tauCaseTypes,tauCaseNames,tauCaseColors);
            %%figSettings = SIM21Analysis.initFigSettings('','','$\nu\ [MHz]$','$\frac{dT_b}{d\nu}\ [\frac{mK}{MHz}]$');
            %%figSettings.legend = false;
            %%obj.plotStackedFills(figSettings,[outputPath,'T21cmDiffStacked.png'],stackCaseTypes,stackCaseNames,stackCaseColors,1);
            %%obj.plotStackedFills(figSettings,[outputPath,'Tau_T21cmDiffStacked.png'],tauCaseTypes,tauCaseNames,tauCaseColors,1);
        end

        function plotFinalSigScatter(obj,xfield,yfield,figSettings,outputName,caseTypes,caseNames,caseColors,caseShapes,caseFilled,freq)
            if ~exist('freq','var')
                freq = 0;
            end
            scatters = {};

            for caseTypeNum = 1:length(caseTypes)
                typeCaseNums = caseTypes{caseTypeNum};
                
                if ~ isempty(typeCaseNums)
                    x = obj.getField(typeCaseNums,xfield);
                    if xfield{end} == 'z'
                        x = x+1;
                    end
                    y = obj.getField(typeCaseNums,yfield);
                    if yfield{end} == 'z'
                        y = y+1;
                    end
                    scatters{end+1} = SIM21Analysis.plotScatter([x;y],caseNames{caseTypeNum},caseColors{caseTypeNum},caseShapes{caseTypeNum},caseFilled{caseTypeNum});
                end
            end

            figSettings.plots = [scatters, figSettings.plots];
            if freq
                for ind = 1:length(figSettings.plots)
                    figSettings.plots{ind}.x = SIM21Utils.T21cmWaveLength./(figSettings.plots{ind}.x);
                end
            end
            if ~isfield(figSettings,'xTick')
                figSettings.xTick = [20:20:400];
            end
            f = figure();
            SIM21Analysis.plotPlot(f,figSettings);
            if freq
                obj.addFrequencyZAxis(1.3.^[0:5]);
            else
                box on;
            end

            saveas(f,outputName);
            saveas(f,outputName(1:end-4),'epsc');
            saveas(f,outputName(1:end-4),'fig');
        end


        function plotStackedFills(obj,figSettings,outputName,caseTypes,caseNames,caseColors,d)
            if ~ exist('d','var')
                d=0;
            end
            for caseTypeInd = 1:length(caseTypes)
                typeCaseNums = caseTypes{caseTypeInd};

                if ~ isempty(typeCaseNums)
                    T21cmDatas = [];
                    for caseNum = typeCaseNums;
                        T21cmData = SIM21Analysis.getZData(obj.paramCases(caseNum).c,'T21cm');
                        T21cmData = [SIM21Utils.T21cmWaveLength./(T21cmData(1,:)+1); T21cmData(2,:)];
                        if d
                            T21cmData = [T21cmData(1,1:end-1);diff(T21cmData(2,:)) ./ diff(T21cmData(1,:))];
                        end
                        T21cmDatas = [T21cmDatas; T21cmData(2,:)];
                    end
                    if length(typeCaseNums) == 1
                        figSettings.plots{end+1} = SIM21Analysis.plotLine([T21cmData(1,:);T21cmData(2,:)],caseNames{caseTypeInd},caseColors{caseTypeInd},'-',2);
                    else
                        figSettings.plots{end+1} = SIM21Analysis.plotFill([T21cmData(1,:); min(T21cmDatas)],[T21cmData(1,:); max(T21cmDatas)],caseNames{caseTypeInd},caseColors{caseTypeInd},0.5);
                    end
                end
            end

            f = figure();
            figSettings.xTick = [0:30:300];
            SIM21Analysis.plotPlot(f,figSettings);
            obj.addFrequencyZAxis(1.5.^[0:8]);

            saveas(f,outputName);
            saveas(f,outputName(1:end-4),'epsc');
        end

        function addFrequencyZAxis(obj,tickVals)
            axnu = gca;
            axz = axes('Position',axnu.Position,'XAxisLocation','top','YAxisLocation','right','Color','none');
            axz.YTick = [];
            axnuXLim = axnu.XLim;
            gca;
            xlim(axnuXLim);
            maxZ1 = SIM21Utils.T21cmWaveLength/axnuXLim(1);
            minZ1 = SIM21Utils.T21cmWaveLength/axnuXLim(2);
            xTickVals = flip(floor(ceil(minZ1).*(tickVals)));
            axz.XTick = SIM21Utils.T21cmWaveLength./xTickVals;
            axz.XTickLabel = num2str(xTickVals');
            xlabel('$1+z$','FontSize',12,'Interpreter','LaTex');
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