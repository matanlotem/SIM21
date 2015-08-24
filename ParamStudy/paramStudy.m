classdef paramStudy < handle
    properties
        dataPath;
        tmpDataPath;
        outputPath;
        workCases;
        specialParams;
        paramCases = struct();
    end
    
    
    methods
        function obj = paramStudy()
            addpath('/a/home/cc/tree/taucc/students/physics/matanlotem/Work/SIM21/lib/');
            
            %paramDataPath = '/a/home/cc/tree/taucc/students/physics/matanlotem/Work/ParamStudy/ParamStudy.txt';
            paramDataPath = 'ParamStudy.txt';
            cubeNum = 9;
            obj.dataPath = '/scratch300/matanlotem/Data/';
            obj.tmpDataPath = '/scratch/matanlotem/Data/';
            obj.outputPath = '/scratch300/matanlotem/ParamStudy/';
            obj.workCases = [1:33,52:65];
            %obj.workCases = [1];
            
            obj.paramCases = obj.getCases(paramDataPath,cubeNum);
            obj.specialParams = obj.getSpecialParams();
            %obj.plotSomething2()
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
        
        
        function specialParams = getSpecialParams(obj)
            for i = obj.workCases
                paramMatName = [obj.getCaseOutputPath(i),'specialParams',obj.paramCases(i).ID,'.mat'];
                
                % Check if specialParams exists and load
                if exist(paramMatName, 'file') == 2
                    specialParams(i)=importdata(paramMatName);
                else
                    disp(i);
                    specialParams(i) = SIM21Analysis.calcSpecialParams(obj.dataPath,obj.getCaseOutputPath(i),obj.paramCases(i).ID);
                end
            end
        end
        
        
        function specialParams = calcSpecialParams(obj)
            % Recalculate specialParams
            for caseNum = obj.workCases
                disp(caseNum);
                specialParams(caseNum) = SIM21Analysis.calcSpecialParams(obj.dataPath,obj.getCaseOutputPath(caseNum),obj.paramCases(caseNum).ID);
            end
        end
        
        
        function plotCaseGraphs(obj)
            for caseNum = obj.workCases
                disp(caseNum);
                SIM21Analysis.plotZGraphs(obj.dataPath,obj.getCaseOutputPath(caseNum),obj.paramCases(caseNum).ID,['Case ',num2str(caseNum)]);
            end
        end
        
        
        function plotEpsGraph(obj,caseNum)
            SIM21Analysis.plotZEpsGraph(obj.tmpDataPath,obj.getCaseOutputPath(caseNum),obj.paramCases(caseNum).ID,['Case ',num2str(caseNum)]);
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
        
        
        function caseOutputPath = getCaseOutputPath(obj,caseNum)
            caseOutputPath = [obj.outputPath,'Case_',num2str(caseNum),'/'];
            [status,message,messageid] = mkdir(caseOutputPath);
        end
    end
end