classdef paramStudy < handle
    properties
        pathExt;
        outputPath;
        regularCase;
        smallVarCases;
        largeVarCases;
        workCases;
        specialParams;
        tempParams;
        paramCases;
    end
    
    
    methods
        function obj = paramStudy()
            addpath('/a/home/cc/tree/taucc/students/physics/matanlotem/Work/SIM21/lib/');
            
            paramDataPath = 'ParamStudy.txt';
            cubeNum = 9;
            obj.pathExt = '';
            obj.outputPath = '/scratch300/matanlotem/ParamStudy/';
            obj.regularCase = [1];
            obj.smallVarCases = [2:33];
            obj.largeVarCases = [52:65];
            obj.largeVarCases = [52:65];
            obj.workCases = [obj.regularCase,obj.smallVarCases,obj.largeVarCases];
            %obj.workCases = [50,51];
            
            obj.paramCases = obj.getCases(paramDataPath,cubeNum);
            obj.specialParams = obj.getSpecialParams();
            obj.tempParams = obj.getTempParams();
            %obj.plotSomething3();
        end


        function paramCases = getCases(obj,paramDataPath,cubeNum)
            rawData = readtable(paramDataPath,'Delimiter','\t');
            %paramCases = struct();
            for i = 1:height(rawData)
                caseNum = rawData.CASE(i);
                paramCases(caseNum).caseNum = caseNum;
                caseName = ['Case ',num2str(caseNum)];
                zeta = char(rawData.Zeta(i));
                if isequal(zeta,'PROBLEM')
                    zeta = NaN;
                    paramCases(caseNum).fxHI = NaN;
                    paramCases(caseNum).atau = NaN;
                    paramCases(caseNum).isgood = false;
                else
                    zeta = str2num(zeta);
                    paramCases(caseNum).fxHI = rawData.xHI_z_6_(i);
                    paramCases(caseNum).atau = rawData.tau_1(i);
                    paramCases(caseNum).isgood = true;
                end
                paramCases(caseNum).c = SIM21Case(caseName,obj.pathExt,obj.getCaseOutputPath(caseNum),...
                                                  cubeNum,rawData.fstar(i),rawData.vbc(i),rawData.vc(i),rawData.fx(i),rawData.sed(i),rawData.tau(i),...
                                                  rawData.LWFlag(i),rawData.LWW_s(i),rawData.pop(i),rawData.func(i),rawData.PH(i),zeta);
            end
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
                xHIData = SIM21Analysis.getZData(obj.paramCases(caseNum).c,SIM21Utils.xHI);
                fxHI = xHIData(2,1);
                xHIData = [];
                atau = SIM21Analysis.checkTau(c);
                
                rowStr = [formatCell(caseNum)                                     ,'\t',...
                          formatCell(obj.paramCases(caseNum).ID)                  ,'\t',...
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
        
        
        function plotAllZGraphs(obj)
            for caseNum = obj.workCases
                disp(caseNum);
                c = obj.paramCases(caseNum).c;
                SIM21Analysis.plotGraphsByZ(c);
                SIM21Analysis.plotEpsByZ(c);
                SIM21Analysis.plotXeByZ(c);
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

        
        function plotSomething(obj)
            f=figure();
            hold on;
            
            minSlope = [obj.specialParams.minSlope];
            scatter([minSlope.z],[minSlope.T],10,'filled');
            
            dx = rand(length(minSlope),1)-0.5;
            dy = rand(length(minSlope),1)-0.5;
            text([minSlope.z]+dx',[minSlope.T]+dy',cellstr(num2str([1:length(minSlope)]')));
            
            hold off;
            outputName = [obj.outputPath,'Graphs/something.png'];
            saveas(f,outputName);
        end
        
        
        function plotSomething2(obj)
            f=figure();
            hold on;
            
            fstarVals = [0.005,0.015,0.05,0.158,0.5];
            vcVals = [4.2,16.5,35.5,76.5];
            fxVals = [0.1,0.16,1,1.58,10];
            sedVals = [1,2,4,5,6];
            tauVals = [0.066,0.082,0.098];
            feedbackFlag = [0,1];
            fsfuncFlag = [1,2];
            
            function drawLine(x,y)
                plot(x,y);
            end
            
            function attributeVal = setAttribute(vals,attrs,val)
                ind = find(vals==val);
                if isempty(ind)
                    % default case
                    attributeVal = attrs(end);
                elseif ind > length(attrs)
                    % default case
                    attributeVal = attrs(end);
                else
                    attributeVal = attrs(ind);
                end
            end
            
            
            for i = obj.workCases
                faceColor = 'b';
                edgeColor = 'k';
                markerShape = 'o';
                markerSize = 40;
                lineWidth = 0.5;
                
                x = obj.specialParams(i).maxSlope.z;
                y = obj.specialParams(i).maxSlope.T;
                
                %faceColor = setAttribute(tauVals,['r','g','y','b'],obj.paramCases(i).tau);
                %faceColor = setAttribute(fstarVals,['b','g','y','r','c','b'],obj.paramCases(i).fstar);
                faceColor = setAttribute(fxVals,['b','g','y','r','c','b'],obj.paramCases(i).fx);
                %fxVals = [0.1,0.16,1,1.58,10];
                
                %markerShape = setAttribute(fxVals,['h','d','o','s','h','o'],obj.paramCases(i).fx);
                %markerShape = setAttribute(vcVals,['s','h','d','p','o'],obj.paramCases(i).vc);
                
                %edgeColor = setAttribute(sedVals,['k','k','m','k'],obj.paramCases(i).sed);
                %markerSize = obj.specialParams(i).maxSlope.slope;
                
                scatter(x,y,markerSize,markerShape,...
                       'MarkerEdgeColor',edgeColor,...
                       'MarkerFaceColor',faceColor,...
                       'LineWidth',lineWidth);
            end
            
            iWorkCases = sum(bsxfun(@eq,[1:65],obj.workCases'));
            for vcVal = vcVals
                ind = find([obj.paramCases.vc] == vcVal .* iWorkCases);
                maxSlope = [obj.specialParams(ind).maxSlope];
                plot([maxSlope.z],[maxSlope.T]);
            end
            
            hold off;
            outputName = [obj.outputPath,'Graphs/something2.png'];
            saveas(f,outputName);
        end
        

        function plotSomething3(obj)
            outputPath = [obj.outputPath,'Graphs/'];
            figSettings.xLabel = 'z';
            figSettings.yLabel = 'T';
            figSettings.title = 'Min T21cm - T(z)';
            obj.plotSigScatter({'specialParams','minT21cm','z'},{'specialParams','minT21cm','T'},figSettings,[outputPath,'minT21cm.png']);
            figSettings.title = 'Max T21cm - T(z)';
            obj.plotSigScatter({'specialParams','maxT21cm','z'},{'specialParams','maxT21cm','T'},figSettings,[outputPath,'maxT21cm.png']);

            figSettings.yLabel = 'Slope';
            figSettings.title = 'Min Slope - slope(z)';
            obj.plotSigScatter({'specialParams','minSlope','z'},{'specialParams','minSlope','slope'},figSettings,[outputPath,'minSlope.png']);
            figSettings.title = 'Max Slope - slope(z)';
            obj.plotSigScatter({'specialParams','maxSlope','z'},{'specialParams','maxSlope','slope'},figSettings,[outputPath,'maxSlope.png']);

            figSettings.xLabel = 'z Min';
            figSettings.yLabel = 'z Max';
            figSettings.title = 'Max T - Min T - zMax(zMin)';
            obj.plotSigScatter({'specialParams','minT21cm','z'},{'specialParams','maxT21cm','z'},figSettings,[outputPath,'minMaxT21cm.png']);
            figSettings.title = 'Max Slope - Min Slope - zMax(zMin)';
            obj.plotSigScatter({'specialParams','minSlope','z'},{'specialParams','maxSlope','z'},figSettings,[outputPath,'minMaxSlope.png']);    

            figSettings.xLabel = 'zeta';
            figSettings.yLabel = 'xHI50 z';
            figSettings.title = 'xHI50 Z(zeta)';
            obj.plotSigScatter({'paramCases','zeta'},{'specialParams','xHI50','z'},figSettings,[outputPath,'xHI50zeta.png']);    
        end


        function plotSomething4(obj)
            outputPath = [obj.outputPath,'Graphs/'];

            xfield = {'specialParams','minSlope','z'};
            yfield = {'specialParams','maxSlope','z'};
            figSettings.title = 'Max Slope - Min Slope - zMax(zMin)';
            %figSettings.lines = [obj.getLines(xfield,yfield,{'paramCases','vc'})];
            %figSettings.lines = [obj.getLines(xfield,yfield,{'paramCases','fx'})];
            figSettings.lines = [obj.getLines(xfield,yfield,{'paramCases','fstar'})];
            obj.plotSigScatter(xfield,yfield,figSettings,[outputPath,'minMaxSlope.png']);

            figSettings.title = 'Max Slope - vc,fstarVals';
            xfield = {'paramCases','vc'};
            yfield = {'specialParams','maxSlope','z'};
            figSettings.lines = [obj.getLines(xfield,yfield,{'paramCases','fstar'})];
            obj.plotSigScatter(xfield,yfield,figSettings,[outputPath,'MaxSlope - vc,fstar.png']);

            figSettings = struct();
            xfield = {'paramCases','vc'};
            yfield = {'tempParams','xHI75_25','zDiff'};
            obj.plotSigScatter(xfield,yfield,figSettings,[outputPath,'xHI75-xHI25 - vc.png']);

            figSettings = struct();
            xfield = {'paramCases','fstar'};
            yfield = {'specialParams','maxSlope','z'};
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

            outputPath = [obj.outputPath,'Graphs/findcorr',num2str(corrLevel),'_'];
            for i = 1:length(row)
                xfield = corrParams{row(i)};
                yfield = corrParams{col(i)};
                xLabel = strjoin(xfield(2:end),'.');
                yLabel = strjoin(yfield(2:end),'.');
                figSettings.xLabel = xLabel;
                figSettings.yLabel = yLabel;
                figSettings.title = [yLabel,' (',xLabel,') - ',num2str(corrs(indCorr(i)))];
                obj.plotSigScatter(xfield,yfield,figSettings,[outputPath,yLabel,'-',xLabel,'.png']);
            end
        end


        function plotSigScatter(obj,xfield,yfield,figSettings,outputName)
            f=figure();
            hold on

            caseTypes = {obj.regularCase,obj.smallVarCases,obj.largeVarCases};
            %caseColors = ['r','g','b'];
            for caseType = 1:length(caseTypes)
                cases = caseTypes{caseType};
                x = obj.getField(cases,xfield);
                y = obj.getField(cases,yfield);
                %scatter(x,y,caseColors(caseType))
                scatter(x,y)
            end
            
            ax = gca;
            if isfield(figSettings,'log')
                if max(figSettings == 'x')
                    ax.XScale = 'log';
                end
                if max(figSettings == 'y')
                    ax.YScale = 'log';
                end
            end

            % Interesting lines (for example zero line)
            if isfield(figSettings,'lines')
                for figLine = figSettings.lines
                    plot(figLine{1}(1,:),figLine{1}(2,:));
                end
            end

            % Titles and Labels
            xLabel = strjoin(xfield(2:end),'.');
            yLabel = strjoin(yfield(2:end),'.');
            if isfield(figSettings,'title');
                title(figSettings.title,'FontSize',18);
            else
                title([yLabel,' (',xLabel,')'],'FontSize',18);
            end
            if isfield(figSettings,'xLabel');
                xlabel(figSettings.xLabel,'FontSize',12);
            else
                xlabel(xLabel);
            end
            if isfield(figSettings,'yLabel');
                ylabel(figSettings.yLabel,'FontSize',12);
            else
                ylabel(yLabel);
            end

                

            legend('regelur','small variation','large variation','Location','bestoutside');

            saveas(f,outputName);
        end


        function lines = getLines(obj,xfield,yfield,linefield)
            lineVals = obj.getField(obj.workCases,linefield);
            uLineVals = unique(lineVals);
            for i = 1:length(uLineVals)
                cases = obj.workCases(lineVals==uLineVals(i));
                lines{i} = cat(1,obj.getField(cases,xfield),obj.getField(cases,yfield));
            end
        end


        function vals = getField(obj,cases,fieldList)
            vals = zeros(1,length(cases));
            for i = 1:length(vals)
                lst = getfield(obj,fieldList{1});
                val = lst(cases(i));

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