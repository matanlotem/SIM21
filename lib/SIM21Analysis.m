classdef SIM21Analysis
    properties(Constant)
        PwSpZ = [6:0.1:15,16:40];
    end
    
    methods(Static)
        function plotAllGraphsByZ(cases)
            SIM21Analysis.plotGraphsByZ(cases);
            SIM21Analysis.plotEpsByZ(cases);
            SIM21Analysis.plotXeByZ(cases);
        end


        function plotGraphsByZ(cases,varargin)
            % Plot TK, T21cm and xHI Graphs
            SIM21Analysis.plotTKByZ(cases,varargin{:});
            SIM21Analysis.plotT21cmByZ(cases,varargin{:});
            SIM21Analysis.plotXHIByZ(cases,varargin{:});
        end
        

        function plotTKByZ(cases,varargin)
            TK = SIM21Utils.dataTypes.TK;
            [figName,figSettings] = SIM21Analysis.initCaseByZFig(cases,TK,'TK(z)','z+1','TK',varargin{:});
            figSettings.log = 'xy';
            figSettings.xLim = [min(TK.z),max(TK.z)] + 1;
            figSettings.yTick = [3,10,30,100,300,1000,3000,10000];
            figSettings.plots{end+1} = SIM21Analysis.plotLine(SIM21Analysis.zPlus1(cat(1,TK.z,SIM21Gets.getTcmb(TK.z))),'TCMB','k',':',0.2);
            SIM21Analysis.plotData(figName,figSettings);
        end


        function plotT21cmByZ(cases,varargin)
            T21cm = SIM21Utils.dataTypes.T21cm;
            [figName,figSettings] = SIM21Analysis.initCaseByZFig(cases,T21cm,'T21cm(z)','z+1','T21cm',varargin{:});
            figSettings.log = 'x';
            figSettings.xLim = [min(T21cm.z),max(T21cm.z)] + 1;
            figSettings.plots{end+1} = SIM21Analysis.plotLine(SIM21Analysis.zPlus1(cat(1,T21cm.z,zeros(1,length(T21cm.z)))),'T0','k',':',0.2);
            SIM21Analysis.plotData(figName,figSettings);
        end


        function plotXHIByZ(cases,varargin)
            xHI = SIM21Utils.dataTypes.xHI;
            [figName,figSettings] = SIM21Analysis.initCaseByZFig(cases,xHI,'xHI(z)','z+1','xHI',varargin{:});
            figSettings.log = 'x';
            figSettings.xLim = [min(xHI.z),max(xHI.z)] + 1;
            figSettings.yLim = [0,1.1];
            SIM21Analysis.plotData(figName,figSettings);
        end
        
        
        function plotEpsByZ(cases,varargin)
            eps = SIM21Utils.dataTypes.eps;
            [figName,figSettings] = SIM21Analysis.initCaseByZFig(cases,eps,'eps(z)','z+1','eps',varargin{:});
            figSettings.log = 'xy';
            figSettings.yTick = repmat([10],1,13).^[-44:3:-8];
            SIM21Analysis.plotData(figName,figSettings);
        end


        function plotXeByZ(cases,varargin)
            xe = SIM21Utils.dataTypes.xe;
            [figName,figSettings] = SIM21Analysis.initCaseByZFig(cases,xe,'xe(z)','z+1','xe',varargin{:});
            figSettings.log = 'x';
            figSettings.xLim = [min(xe.z),max(xe.z)] + 1;
            figSettings.yLim = [0,1.1];
            SIM21Analysis.plotData(figName,figSettings);
        end

        function plotXeXHIByZ(c)
            xe = SIM21Utils.dataTypes.xe;
            xHI = SIM21Utils.dataTypes.xHI;
            figName = [c.outputPath,'xHI-xe',c.ID,'.png'];
            figSettings = SIM21Analysis.initFigSettings('xHI & 1-xe',c.name,'z+1','x');
            figSettings.log = 'x';
            figSettings.plots{end+1} = SIM21Analysis.plotLine(SIM21Analysis.zPlus1(SIM21Analysis.getZData(c,xHI)),'xHI');
            xeData = SIM21Analysis.getZData(c,xe);
            xe1Data = cat(1,xeData(1,:),1-xeData(2,:));
            figSettings.plots{end+1} = SIM21Analysis.plotLine(SIM21Analysis.zPlus1(xe1Data),'1-xe');

            figSettings.xLim = [min(xHI.z),max(xHI.z)] + 1;
            figSettings.yLim = [0,1.1];

            SIM21Analysis.plotData(figName,figSettings);
        end

        function plotTKTSByZ(cases,varargin)
            TK = SIM21Utils.dataTypes.TK;
            TS = SIM21Utils.dataTypes.TS;
            [figName,figSettings] = SIM21Analysis.initCaseByZFig(cases,TS,'T(z)','z+1','TS',varargin{:});
            figSettings.log = 'xy';
            figSettings.xLim = [min(TK.z),max(TK.z)] + 1;
            figSettings.yTick = [3,10,30,100,300,1000,3000,10000];
            figSettings.plots{end+1} = SIM21Analysis.plotLine(SIM21Analysis.zPlus1(SIM21Analysis.getZData(cases(1),TK)),'TK');
            figSettings.plots{end+1} = SIM21Analysis.plotLine(SIM21Analysis.zPlus1(cat(1,TK.z,SIM21Gets.getTcmb(TK.z))),'TCMB','k',':',0.2);
            SIM21Analysis.plotData(figName,figSettings);
        end


        
        function plotSanityGraphs(c,varargin)
            f = figure();
            xHI = SIM21Utils.dataTypes.xHI;
            xe = SIM21Utils.dataTypes.xe;
            T21cm = SIM21Utils.dataTypes.T21cm;
            TK = SIM21Utils.dataTypes.TK;
            eps = SIM21Utils.dataTypes.eps;
            
            xHIData = SIM21Analysis.getZData(c,'xHI');
            reionZ = xHIData(1,max(find(xHIData(2,:)<0.06)));

            % T21cm
            [figName,figSettings] = SIM21Analysis.initCaseByZFig(c,T21cm,'T21cm(z)','z+1','T21cm',varargin{:});
            figSettings.log = 'x';
            figSettings.plots{end+1} = SIM21Analysis.plotLine(SIM21Analysis.zPlus1(cat(1,T21cm.z,zeros(1,length(T21cm.z)))),'T0','k',':',0.2);

            figSettings.xLim = [min(T21cm.z),max(T21cm.z)] + 1;
            figSettings.plots{end+1} = SIM21Analysis.plotLine([reionZ+1,reionZ+1;min(figSettings.plots{1}.y),max(figSettings.plots{1}.y)],'Reion','k','--',0.2);
            subplot(2,2,1)
            SIM21Analysis.plotPlot(f,figSettings);
            title(figSettings.title,'FontSize',12);
            xlabel(figSettings.xLabel,'FontSize',9);
            ylabel(figSettings.yLabel,'FontSize',9);
            legend([]);

            % TK
            [figName,figSettings] = SIM21Analysis.initCaseByZFig(c,TK,'TK(z)','z+1','TK',varargin{:});
            figSettings.log = 'xy';
            figSettings.xLim = [min(TK.z),max(TK.z)] + 1;
            figSettings.yTick = [3,10,30,100,300,1000,3000,10000];
            figSettings.plots{end+1} = SIM21Analysis.plotLine(SIM21Analysis.zPlus1(cat(1,TK.z,SIM21Gets.getTcmb(TK.z))),'TCMB','k',':',0.2);

            figSettings.xLim = [min(T21cm.z),max(T21cm.z)] + 1;
            figSettings.plots{end+1} = SIM21Analysis.plotLine([reionZ+1,reionZ+1;min(figSettings.plots{1}.y),max(figSettings.plots{1}.y)],'Reion','k','--',0.2);
            subplot(2,2,2)
            SIM21Analysis.plotPlot(f,figSettings);
            title(figSettings.title,'FontSize',12);
            xlabel(figSettings.xLabel,'FontSize',9);
            ylabel(figSettings.yLabel,'FontSize',9);
            legend([]);
            
            % xHI-xe
            figSettings = SIM21Analysis.initFigSettings('xHI & 1-xe (z)','','z+1','x');
            figSettings.plots{end+1} = SIM21Analysis.plotLine(SIM21Analysis.zPlus1(SIM21Analysis.getZData(c,xHI)),'xHI');
            xeData = SIM21Analysis.getZData(c,xe);
            xe1Data = cat(1,xeData(1,:),1-xeData(2,:));
            figSettings.plots{end+1} = SIM21Analysis.plotLine(SIM21Analysis.zPlus1(xe1Data),'1-xe');
            
            figSettings.log = 'x';
            figSettings.xLim = [min(xHI.z),max(xHI.z)] + 1;
            figSettings.yLim = [0,1.1];

            figSettings.xLim = [min(T21cm.z),max(T21cm.z)] + 1;
            figSettings.plots{end+1} = SIM21Analysis.plotLine([reionZ+1,reionZ+1;min(figSettings.plots{1}.y),max(figSettings.plots{1}.y)],'Reion','k','--',0.2);
            subplot(2,2,3)
            SIM21Analysis.plotPlot(f,figSettings);
            title(figSettings.title,'FontSize',12);
            xlabel(figSettings.xLabel,'FontSize',9);
            ylabel(figSettings.yLabel,'FontSize',9);
            legend([]);

            % eps
            [figName,figSettings] = SIM21Analysis.initCaseByZFig(c,eps,'eps(z)','z+1','eps',varargin{:});
            figSettings.log = 'xy';
            figSettings.yTick = repmat([10],1,13).^[-44:3:-8];

            figSettings.xLim = [min(T21cm.z),max(T21cm.z)] + 1;
            figSettings.plots{end+1} = SIM21Analysis.plotLine([reionZ+1,reionZ+1;min(figSettings.plots{1}.y),max(figSettings.plots{1}.y)],'Reion','k','--',0.2);
            subplot(2,2,4)
            SIM21Analysis.plotPlot(f,figSettings);
            title(figSettings.title,'FontSize',12);
            xlabel(figSettings.xLabel,'FontSize',9);
            ylabel(figSettings.yLabel,'FontSize',9);
            legend([]);

            outputName = [c.outputPath,'Sanity-',c.name(6:end),c.ID,'.png'];
            saveas(f,outputName);
        end


        function [figName,figSettings] = initCaseByZFig(cases,dataType,title,xLabel,yLabel,varargin)
            c = cases(end);
            if length(varargin)==1
                figName = varargin{1};
                add2Title = '';
            elseif length(varargin) > 1
                figName = [varargin{1},dataType.magic,'_',varargin{2},'.png'];
                add2Title = varargin{2};
            else
                figName = [c.outputPath,dataType.magic,c.ID,'.png'];
                add2Title = c.name;
            end
            figSettings = SIM21Analysis.initFigSettings(title,add2Title,xLabel,yLabel);
            for c = cases
                figSettings.plots{end+1} = SIM21Analysis.plotLine(SIM21Analysis.zPlus1(SIM21Analysis.getZData(c,dataType)),c.name);
            end
        end


        function figSettings = initFigSettings(title,add2Title,xLabel,yLabel)
            if ~isempty(add2Title)
                add2Title = [' - ',add2Title];
            end
            figSettings.title = [title,add2Title];
            figSettings.xLabel = xLabel;
            figSettings.yLabel = yLabel;
            figSettings.plots = {};
            figSettings.legend = true;
        end
        

        function dataMat = zPlus1(dataMat)
            dataMat(1,:) = dataMat(1,:) + 1;
        end
        
        function dataMat = getZData(c,dataType)
            % Create or import mean data matrix
            dataType = SIM21Utils.getDataType(dataType);
            dataName = [c.outputPath,dataType.magic,'_Data',c.ID,'.mat'];
            
            % Check if output exists and load
            if exist(dataName, 'file') == 2
                dataMat=importdata(dataName);
            else
                % Calculate Mean
                SIM21Analysis.message(['calculating ',dataType.magic,' mean']);

                if isequal(dataType,SIM21Utils.getDataType('TS'));
                    dataMat = SIM21Analysis.getTsZData(c);
                else
                    dataMat=cat(1,dataType.z,nanmean(nanmean(nanmean(c.getData(dataType),4),3),2)');
                end
                
                % Save Output
                SIM21Analysis.message('saving data');
                save(dataName,'dataMat');
            end
        end


        function dataMat = getTsZData(c)
            TKData = c.getData('TK');
            TS = SIM21Utils.getDataType('TS');
            
            global ID
            global pathname_Data1
            global pathname_DataBackgrounds
            global delta_cube

            pathname_Data1 = c.tmpDataPath;
            pathname_DataBackgrounds = SIM21Utils.paths.dataBackgrounds;
            ID = c.ID;
            pathname_Data1 = c.tmpDataPath;
            delta_cube=importdata(strcat(pathname_DataBackgrounds,'my',num2str(c.ncube),'_d.dat'));

            TSMean = [];
            for z = TS.z
                TSData = SIM21Gets.getTs(squeeze(TKData(find(SIM21Utils.dataTypes.TK.z == z),:,:,:)),z,c.ncube,c.fstar,c.vbc,c.vc,c.fx,c.sed,c.zeta,c.feedback,c.delayParam,c.pop,c.fsfunc,c.phVersion);
                save(SIM21Utils.getDataFileName(c,'TS',z),'TSData');
                TSMean = [TSMean,[nanmean(nanmean(nanmean(TSData)))]];
            end
            dataMat=cat(1,TS.z,TSMean);
        end


        function XYData=interpData(dataMat,interpStep)
            % Interpolate data matrix
            SIM21Analysis.message('interpolating');
            x = dataMat(1,:);
            y = dataMat(2,:);
            x_i = min(x):interpStep:max(x);
            y_i = interp1(x,y,x_i,'spline');
            XYData = cat(1,x_i,y_i);
        end
        


        function plotData(outputName,figSettings)
            % Plot data
            f=figure();
            SIM21Analysis.plotPlot(f,figSettings);
            saveas(f,outputName);
        end

        function f = plotPlot(f,figSettings)
            hold on;
            
            plotNames = {};
            for figPlot = figSettings.plots
                figPlot = figPlot{1};
                switch figPlot.type
                case 'line'
                    h = plot(figPlot.x,figPlot.y);
                    if isfield(figPlot,'lineColor')
                        h.Color = figPlot.lineColor;
                    end
                    if isfield(figPlot,'lineStyle')
                        h.LineStyle = figPlot.lineStyle;
                    end
                    if isfield(figPlot,'lineWidth')
                        h.LineWidth = figPlot.lineWidth;
                    end
                case 'scatter'
                    h = scatter(figPlot.x,figPlot.y);
                    if isfield(figPlot,'color')
                        h.MarkerEdgeColor = figPlot.color;
                    else
                        h.MarkerEdgeColor = h.CData;
                    end
                    if isfield(figPlot,'shape')
                        h.Marker = figPlot.shape;
                    end
                case 'fill'
                    if isfield(figPlot,'color')
                        fillColor = figPlot.color;
                    else
                        fillColor = 'k';
                    end
                    if isfield(figPlot,'alpha')
                        fillAlpha = figPlot.alpha;
                    else
                        fillAlpha = 1;
                    end
                    h = fill(figPlot.x,figPlot.y,fillColor,'EdgeColor','None','facealpha',fillAlpha);
                end
                plotNames{end+1} = figPlot.name;
            end

            ax = gca;

            % Set logarithmic axes
            if isfield(figSettings,'log')
                if max(figSettings.log == 'x')
                    ax.XScale = 'log';
                end
                if max(figSettings.log == 'y')
                    ax.YScale = 'log';
                end
            end

            % Set limits
            if isfield(figSettings,'xLim')
                xlim(figSettings.xLim);
            end
            if isfield(figSettings,'yLim')
                ylim(figSettings.yLim);
            end

            % Set tick marks
            if isfield(figSettings,'xTick');
                ax.XTick = figSettings.xTick;
            else
                1+1;
                ax.XTick = floor(min(ax.XLim)/10)*10:10:max(ax.XLim);
            end
            if isfield(figSettings,'yTick')
                ax.YTick = figSettings.yTick;
            else
                2+2;
                numOfYTicks = 5;
                yhop = (max(ax.YLim)-min(ax.YLim))/numOfYTicks;
                ymag = 10^floor(log10(yhop));
                yhop = floor(yhop/ymag)*ymag;
                ax.YTick = (ceil(min(ax.YLim) / yhop) * yhop):yhop:max(ax.YLim);
            end
            
            % Titles and Labels
            if isfield(figSettings,'title');
                title(figSettings.title,'FontSize',18,'Interpreter','LaTex');
            end
            if isfield(figSettings,'xLabel');
                xlabel(figSettings.xLabel,'FontSize',12,'Interpreter','LaTex');
            end
            if isfield(figSettings,'yLabel');
                ylabel(figSettings.yLabel,'FontSize',12,'Interpreter','LaTex');
            end
            hold off;
            
            % Legend
            if figSettings.legend
                legend(plotNames,'Location','bestoutside');
            end
        end


        function pline = plotLine(XYData,name,lineColor,lineStyle,lineWidth)
            % create line object
            pline.type = 'line';
            pline.x = XYData(1,:);
            pline.y = XYData(2,:);
            pline.name = name;
            if exist('lineColor','var')
                pline.lineColor = lineColor;
            end
            if exist('lineStyle','var')
                pline.lineStyle = lineStyle;
            end
            if exist('lineWidth','var')
                pline.lineWidth = lineWidth;
            end
        end

        function pfill = plotFill(XYData1,XYData2,name,Color,Alpha)
            % create line object
            pfill.type = 'fill';
            pfill.x = [XYData1(1,:), flip(XYData2(1,:))];
            pfill.y = [XYData1(2,:), flip(XYData2(2,:))];
            pfill.name = name;
            if exist('Color','var')
                pfill.color = Color;
            end
            if exist('Alpha','var')
                pfill.alpha = Alpha;
            end
        end


        function pscatter = plotScatter(XYData,name,color,shape)
            % create line object
            pscatter.type = 'scatter';
            pscatter.x = XYData(1,:);
            pscatter.y = XYData(2,:);
            pscatter.name = name;
            if exist('color','var')
                pscatter.color = color;
            end
            if exist('shape','var')
                pscatter.shape = shape;
            end
        end

        
        function specialParams=calcSpecialParams(c)
            %Calculate interesting parameters for specific run
            SIM21Analysis.message('calculating parameters');
            interpStep=0.001;
            zMax = 70;
            
            xHIData = SIM21Analysis.interpData(SIM21Analysis.getZData(c,SIM21Utils.dataTypes.xHI),interpStep);
            xHIZMaxInd = find(xHIData(1,:)==min(max(xHIData(1,:)),zMax));
            TKData = SIM21Analysis.interpData(SIM21Analysis.getZData(c,SIM21Utils.dataTypes.TK),interpStep);
            TKZMaxInd = find(TKData(1,:)==min(max(TKData(1,:)),zMax));
            TSData = SIM21Analysis.interpData(SIM21Analysis.getZData(c,SIM21Utils.dataTypes.TS),interpStep);
            TSZMaxInd = find(TSData(1,:)==min(max(TSData(1,:)),zMax));
            T21cmData = SIM21Analysis.interpData(SIM21Analysis.getZData(c,SIM21Utils.dataTypes.T21cm),interpStep);
            T21cmZMaxInd = max(find(diff(T21cmData(2,:))>0))+1;

            
            % MIN / MAX T21cm
            minT21cmInd = find(T21cmData(2,:)==min(T21cmData(2,1:T21cmZMaxInd)));
            maxT21cmInd = find(T21cmData(2,:)==max(T21cmData(2,1:T21cmZMaxInd)));
            specialParams.minT21cm.z = T21cmData(1,minT21cmInd);
            specialParams.minT21cm.T = T21cmData(2,minT21cmInd);
            specialParams.maxT21cm.z = T21cmData(1,maxT21cmInd);
            specialParams.maxT21cm.T = T21cmData(2,maxT21cmInd);
            specialParams.fmaxT21cm.z = T21cmData(1,T21cmZMaxInd);
            specialParams.fmaxT21cm.T = T21cmData(2,T21cmZMaxInd);
            
            % MIN / MAX Slope
            %slope = diff(T21cmData(2,1:T21cmZMaxInd)) ./ diff(T21cmData(1,1:T21cmZMaxInd));
            %minSlopeInd = find(slope==min(slope));
            %maxSlopeInd = find(slope==max(slope));
            %specialParams.minSlope.z = T21cmData(1,minSlopeInd);
            %specialParams.minSlope.T = T21cmData(2,minSlopeInd);
            %specialParams.minSlope.slope = slope(minSlopeInd);
            %specialParams.maxSlope.z = T21cmData(1,maxSlopeInd);
            %specialParams.maxSlope.T = T21cmData(2,maxSlopeInd);
            %specialParams.maxSlope.slope = slope(maxSlopeInd);

            % MIN / MAX Nu Slope
            slope = diff(T21cmData(2,1:T21cmZMaxInd)) ./ diff(SIM21Utils.T21cmWaveLength./(T21cmData(1,1:T21cmZMaxInd)+1));
            minSlopeInd = find(slope==min(slope));
            maxSlopeInd = find(slope==max(slope));
            specialParams.minSlope.z = T21cmData(1,minSlopeInd);
            specialParams.minSlope.nu = SIM21Utils.T21cmWaveLength/(specialParams.minSlope.z+1);
            specialParams.minSlope.T = T21cmData(2,minSlopeInd);
            specialParams.minSlope.slope = slope(minSlopeInd);
            specialParams.maxSlope.z = T21cmData(1,maxSlopeInd);
            specialParams.maxSlope.nu = SIM21Utils.T21cmWaveLength/(specialParams.maxSlope.z+1);
            specialParams.maxSlope.T = T21cmData(2,maxSlopeInd);
            specialParams.maxSlope.slope = slope(maxSlopeInd);
            
            % X Crossing
            xCrossInd = find(diff(sign(T21cmData(2,1:T21cmZMaxInd))));
            specialParams.xCross.z = T21cmData(1,xCrossInd);
            specialParams.xCross.z1 = max(specialParams.xCross.z); % in case there are several zero crossings

            % MIN TK
            minTKInd = find(TKData(2,:)==min(TKData(2,1:TKZMaxInd)));
            specialParams.minTK.z = TKData(1,minTKInd);
            specialParams.minTK.T = TKData(2,minTKInd);

            % MIN TS
            minTSInd = find(TSData(2,:)==min(TSData(2,1:TSZMaxInd)));
            specialParams.minTS.z = TSData(1,minTSInd);
            specialParams.minTS.T = TSData(2,minTSInd);
            specialParams.minTS.xHI = interp1(xHIData(1,:),xHIData(2,:),specialParams.minTS.z,'spline');
            
            % xHI Percentage
            xHI95Ind = find(diff(sign(xHIData(2,1:xHIZMaxInd)-0.95)));
            specialParams.xHI95.z = xHIData(1,xHI95Ind);
            xHI80Ind = find(diff(sign(xHIData(2,1:xHIZMaxInd)-0.80)));
            specialParams.xHI80.z = xHIData(1,xHI80Ind);
            xHI75Ind = find(diff(sign(xHIData(2,1:xHIZMaxInd)-0.75)));
            specialParams.xHI75.z = xHIData(1,xHI75Ind);
            xHI50Ind = find(diff(sign(xHIData(2,1:xHIZMaxInd)-0.50)));
            specialParams.xHI50.z = xHIData(1,xHI50Ind);
            xHI25Ind = find(diff(sign(xHIData(2,1:xHIZMaxInd)-0.25)));
            specialParams.xHI25.z = xHIData(1,xHI25Ind);
            xHI005Ind = find(diff(sign(xHIData(2,1:xHIZMaxInd)-0.05)));
            specialParams.xHI005.z = xHIData(1,xHI005Ind);
            xHI001Ind = find(diff(sign(xHIData(2,1:xHIZMaxInd)-0.01)));
            specialParams.xHI001.z = xHIData(1,xHI001Ind);
            xHI0Ind = max(find(xHIData(2,1:xHIZMaxInd)==0));
            specialParams.xHI0.z = xHIData(1,xHI0Ind);
            
            % Heating Transition: Tcmb = Tk
            TCMBData = SIM21Analysis.interpData(cat(1,SIM21Utils.dataTypes.TK.z,SIM21Gets.getTcmb(SIM21Utils.dataTypes.TK.z)),interpStep);
            THTInd = find(diff(sign(TKData(2,:)-TCMBData(2,:))));
            specialParams.THT.z = TKData(1,THTInd);
            specialParams.THT.T = TKData(2,THTInd);
            % Heating Transition: Tb = 0
            T21cm0Ind = max(find(diff(sign(T21cmData(2,maxT21cmInd:minT21cmInd)))) + maxT21cmInd);
            specialParams.T21cm0.z = T21cmData(1,T21cm0Ind);
            
            % Save matrix
            fileName = [c.outputPath,'specialParams',c.ID];
            save([fileName,'.mat'],'specialParams');
            
            % Save readable output
            fid = fopen([fileName,'.txt'],'w');
            fwrite(fid,SIM21Analysis.niceSpecialParams(specialParams));
            fclose(fid);
        end
        
        
        function niceOutput = niceSpecialParams(specialParams)
            niceOutput = ['Min T\n\tz = ',num2str(specialParams.minT21cm.z),'\n\tT = ',num2str(specialParams.minT21cm.T),'\n',...
                          'Max T\n\tz = ',num2str(specialParams.maxT21cm.z),'\n\tT = ',num2str(specialParams.maxT21cm.T),'\n',... 
                          'Min Slope\n\tz = ',num2str(specialParams.minSlope.z),'\n\tT = ',num2str(specialParams.minSlope.T),'\n\tSlope = ',num2str(specialParams.minSlope.slope),'\n',... 
                          'Max Slope\n\tz = ',num2str(specialParams.maxSlope.z),'\n\tT = ',num2str(specialParams.maxSlope.T),'\n\tSlope = ',num2str(specialParams.maxSlope.slope),'\n',...
                          'Min TK\n\tz = ',num2str(specialParams.minTK.z),'\n\tT = ',num2str(specialParams.minTK.T),'\n',...
                          '0 Crossing\n\tz = ',num2str(specialParams.xCross.z),'\n',...
                          'xHI Percentage\n\t75%% z = ',num2str(specialParams.xHI75.z),'\n\t50%% z = ',num2str(specialParams.xHI50.z),...
                                        '\n\t25%% z = ',num2str(specialParams.xHI25.z),'\n\t0%% z = ',num2str(specialParams.xHI0.z),'\n',...
                          'Heating Transition\n\tz = ',num2str(specialParams.THT.z),'\n\tT = ',num2str(specialParams.THT.T),'\n'];
            niceOutput = sprintf(niceOutput);
        end
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% POWER SPECTRUM STUFF %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        function plotPowerSpectrums(outputPath,runID)
            SIM21Analysis.plotPowerSpectrumByK(outputPath,runID,0.1,['Atomic, Old spectrum - K=',num2str(0.1)],'1+z','k^3P(k)/2\pi^2 [mK^2]');
            SIM21Analysis.plotPowerSpectrumByK(outputPath,runID,0.5,['Atomic, Old spectrum - K=',num2str(0.5)],'1+z','k^3P(k)/2\pi^2 [mK^2]');
            SIM21Analysis.plotPowerSpectrumByZ(outputPath,runID,'Atomic, Old spectrum','k [Mpc^{-1}]','k^3P(k)/2\pi^2 [mK^2]');
        end
        
        
        function plotPowerSpectrumByZ(outputPath,runID,figTitle,figXLabel,figYLabel)
            zs = [8,8.7,10.37,17,11.53,22,30];
            lineColors = {'r','g','b','m','c','k','k'};
            lineStyles = {'-','-','-','-','-','-','--'};
            lineWidths = ones(1,7);
            SIM21Analysis.plotPowerSpectrumByZs(outputPath,runID,zs,lineColors,lineStyles,lineWidths,figTitle,figXLabel,figYLabel);
        end
        
        
        function plotPowerSpectrumByK(outputPath,runID,k,figTitle,figXLabel,figYLabel)
            SIM21Analysis.message('plotting power spectrum');
            
            % Get K data
            MK = SIM21Utils.importMatrix('K');
            % Find index of closest K
            KInd = find(diff(sign(MK-k)));
            
            % Interpolation Parameters
            interpStep = 0.01;
            zPS = min(SIM21Analysis.PwSpZ):interpStep:max(SIM21Analysis.PwSpZ); 
            
            % Plotting Paramters
            figName = [outputPath,'PwSpK_',num2str(k),runID,'.png'];
            delx = 0.04;
            x2=min(SIM21Analysis.PwSpZ):delx:max(SIM21Analysis.PwSpZ); 
            Ind = floor(1+delx*(0:length(x2)-1)/0.01);
            
            f=figure();
            
            function doPowerSpectrum(MName,Msign,lineColor,lineStyle,lineWidth)
                % Import -> Calc for specific K -> Interpolate -> Plot
                PowerMat = importdata([outputPath,MName,runID,'.mat']);
                PS = (MK(KInd).^3.*squeeze(real(PowerMat(:,KInd)))/(2*pi^2));
                PS1 = interp1(SIM21Analysis.PwSpZ,Msign*PS',zPS,'spline');
                loglog(1+x2,PS1(Ind).*(x2>6.9),'Color',lineColor,'LineStyle',lineStyle,'LineWidth',lineWidth);
                hold on
            end
            
            doPowerSpectrum('PowerMat',1,'r','-',1);
            doPowerSpectrum('PowerMat_iso',1,[0,127/255,0],'-',1);
            doPowerSpectrum('PowerMat_X',1,'b','-',1);
            doPowerSpectrum('PowerMat_X',-1,'b',':',2);
            doPowerSpectrum('PowerMat_del',1,'k','-',1);
            
            % Style Plot
            xlim([7,40]);
            ylim([1e-4,1e3]);
            
            title(figTitle,'FontSize',18);
            xlabel(figXLabel,'FontSize',12);
            ylabel(figYLabel,'FontSize',12);
            
            hold off
            
            SIM21Analysis.message('saving plot');
            saveas(f,figName);
        end
        
        
        function plotPowerSpectrumByZs(outputPath,runID,zs,lineColors,lineStyles,lineWidths,figTitle,figXLabel,figYLabel)
            
            MK = SIM21Utils.importMatrix('K');
            PowerMat = importdata([outputPath,'PowerMat',runID,'.mat']);
            
            figName = [outputPath,'PwSpZ',runID,'.png'];
            f=figure();
            
            function doPowerSpectrum(z,lineColor,lineStyle,lineWidth)
                % Calc for specific z -> Interpolate -> Plot
                zInd = find(~(SIM21Analysis.PwSpZ-z));
                if isempty(zInd)
                    zInd = find(diff(sign(SIM21Analysis.PwSpZ-z)));
                end
                
                PS1 = smooth((MK.^3.*squeeze(PowerMat(zInd,:))/(2*pi^2)),'moving');
                loglog(MK,PS1,'Color',lineColor,'LineStyle',lineStyle,'LineWidth',lineWidth);
            end
            
            for i = 1:length(zs)
                doPowerSpectrum(zs(i),lineColors{i},lineStyles{i},lineWidths(i));
                hold on
            end
            
            xlim([0.03,1]);
            ylim([0.01,1000]);
            
            title(figTitle,'FontSize',18);
            xlabel(figXLabel,'FontSize',12);
            ylabel(figYLabel,'FontSize',12);
            
            hold off
            
            SIM21Analysis.message('saving plot');
            saveas(f,figName);
        end
        
        
        %%%%%%%%%%%%%%%%%%%
        %%% OTHER STUFF %%%
        %%%%%%%%%%%%%%%%%%%
        function tauCMB = checkTau(runCase)
            % Aviad's code with fixes - I have no idea what this does...
            load(SIM21Utils.getMatrixPath('Planck_parameters'));

            zion=[0:0.1:15 16:1:20];

            sigmaT = 6.65e-25; % cm^2
            G = 6.67e-8;% cm^3/g s^2
            rcr = 3*(H0*1e5/3e24)^2/(8*pi*G);% gr/cm^3 
            Y = 0.248;
            y = (Y/4)/(1-Y);
            mH = 1.67e-24;% gr
            dldz = c./SIM21Gets.getHz(zion)./(1+zion)*3e24;%cm
            ne = Ob*(1-Y)*(1+y)*(rcr/mH)*(1+zion).^3;%nb./(3e24)^3.*(1+zc2).^3.*(1-xHImatR(indp,:));% 

            xHImat = SIM21Analysis.getZData(runCase,SIM21Utils.dataTypes.xHI);
            xHImat = [zeros(1,60),xHImat(2,:)];
            tauCMB = sum((1-xHImat(1:151)).*ne(1:151).*sigmaT.*dldz(1:151))*0.1+sum((1-xHImat(152:156)).*ne(152:156).*sigmaT.*dldz(152:156))*1;
        end


        function finXHI = checkFinXHI(c)
            xHIData = SIM21Analysis.getZData(c,SIM21Utils.dataTypes.xHI);
            finXHI = xHIData(2,1);
        end
                
        
        function message(msg)
            % SIM21Analysis class output
            if true
                disp(msg);
            end
        end
        
        
        function errorMessage(msg)
            % SIM21Analysis class output
            if true
                disp(msg);
            end
        end
    end
end
