classdef SIM21Analysis
    properties(Constant)
        %Data Matrix settings
        xHIMagic = 'xHI';
        xHIZ = [6:0.1:15,16:60];
        TKMagic = 'TK';
        TKZ = [5:64];
        T21cmMagic = 'T21cm';
        T21cmZ = [6:60];
        NeutMagic = 'Neut';
        NeutZ = [6:0.1:15,16:60];
        EpsMagic = 'eps';
        EpsZ = [6:60];
        XeMagic = 'xe';
        XeZ = [5:64];
        JLWMagic = 'JLW';
        JLWZ = [6:66];
        JalphaMagic = 'Jalpha';
        JalphaZ = [6:60];
        LionMagic = 'Lion';
        LionZ = [6:60];
        PwSpZ = [6:0.1:15,16:40];
    end
    
    methods(Static)
        function analyse()
            dataPath = '/scratch300/matanlotem/Data/';
            outputPath = '/scratch300/matanlotem/Analysis/';
            runID = '_0_0.05_1_16.5_1_1_0.075_0_0_2_1_2';
            
            SIM21Analysis.plotZGraphs(dataPath,outputPath,runID);
            specialParams=SIM21Analysis.calcSpecialParams(dataPath,outputPath,runID);
            SIM21Analysis.plotPowerSpectrums(outputPath,runID);
        end
        
        
        function plotZGraphs(dataPath,outputPath,runID,add2Title)
            % Plot TK, T21cm and xHI Graphs
            SIM21Analysis.plotZTKGraph(dataPath,outputPath,runID,add2Title);
            SIM21Analysis.plotZT21cmGraph(dataPath,outputPath,runID,add2Title);
            SIM21Analysis.plotZXHIGraph(dataPath,outputPath,runID,add2Title);
        end
        
        
        function plotZTKGraph(dataPath,outputPath,runID,add2Title)
            figSettings.log = 'xy';
            figSettings.xLim = [min(SIM21Analysis.TKZ),max(SIM21Analysis.TKZ)+10];
            %figSettings.yLim = [4,1000];
            figSettings.yTick = [3,10,30,100,300,1000,3000,10000];
            if ~isempty(add2Title)
                add2Title = [' - ',add2Title];
            end
            figSettings.title = ['TK(z)',add2Title];
            figSettings.xLabel = '1+z';
            figSettings.yLabel = 'TK [K]';
            figSettings.lines = {cat(1,SIM21Analysis.TKZ(2:end),SIM21Gets.getTcmb(SIM21Analysis.TKZ(2:end)))};
            SIM21Analysis.plotByZ(dataPath,outputPath,SIM21Analysis.TKZ,SIM21Analysis.TKMagic,0,runID,figSettings);
        end
        
        
        function plotZT21cmGraph(dataPath,outputPath,runID,add2Title)
            figSettings.log = 'x';
            figSettings.xLim = [min(SIM21Analysis.T21cmZ),max(SIM21Analysis.T21cmZ)+10];
            if ~isempty(add2Title)
                add2Title = [' - ',add2Title];
            end
            figSettings.title = ['T21cm(z)',add2Title];
            figSettings.xLabel = '1+z';
            figSettings.yLabel = 'T21cm [mK]';
            figSettings.lines = {cat(1,SIM21Analysis.T21cmZ(2:end),zeros(1,length(SIM21Analysis.T21cmZ(2:end))))};
            SIM21Analysis.plotByZ(dataPath,outputPath,SIM21Analysis.T21cmZ,SIM21Analysis.T21cmMagic,0,runID,figSettings);
        end
        
        
        function plotZXHIGraph(dataPath,outputPath,runID,add2Title)
            figSettings.xLim = [min(SIM21Analysis.xHIZ),max(SIM21Analysis.xHIZ)+1];
            figSettings.yLim = [0,1.1];
            if ~isempty(add2Title)
                add2Title = [' - ',add2Title];
            end
            figSettings.title = ['xHI(z)',add2Title];
            figSettings.xLabel = '1+z';
            figSettings.yLabel = 'xHI';
            SIM21Analysis.plotByZ(dataPath,outputPath,SIM21Analysis.xHIZ,SIM21Analysis.xHIMagic,0,runID,figSettings);
        end
        
        
        function plotZEpsGraph(dataPath,outputPath,runID,add2Title)
            figSettings.log = 'y';
            %figSettings.xLim = [min(SIM21Analysis.EpsZ),max(SIM21Analysis.EpsZ)+1];
            figSettings.yTick = repmat([10],1,13).^[-44:3:-8];
            if ~isempty(add2Title)
                add2Title = [' - ',add2Title];
            end
            figSettings.title = ['Epsilon(z)',add2Title];
            figSettings.xLabel = '1+z';
            figSettings.yLabel = 'Epsilon';
            SIM21Analysis.plotByZ(dataPath,outputPath,SIM21Analysis.EpsZ,SIM21Analysis.EpsMagic,0,runID,figSettings);
        end
        
        
        function plotZXeGraph(dataPath,outputPath,runID,add2Title)
            figSettings.log = '';
            %figSettings.yTick = repmat([10],1,13).^[-44:3:-8];
            if ~isempty(add2Title)
                add2Title = [' - ',add2Title];
            end
            figSettings.title = ['XE(z)',add2Title];
            figSettings.xLabel = '1+z';
            figSettings.yLabel = 'XE';
            SIM21Analysis.plotByZ(dataPath,outputPath,SIM21Analysis.XeZ,SIM21Analysis.XeMagic,0,runID,figSettings);
        end
        
        
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
        
        
        function plotByZ(dataPath,outputPath,z,Magic,interpStep,runID,figSettings)
            % Plot different graphs as function of Z
            
            SIM21Analysis.message(['=== STARTING ',Magic,' ===']);
            dataMat = SIM21Analysis.getZData(dataPath,outputPath,z,Magic,runID);
            if interpStep == 0
                XYData = dataMat;
            else
                XYData = SIM21Analysis.interpData(dataMat,interpStep);
            end
            XYData(1,:) = XYData(1,:) + 1; % z+1
            figName = [outputPath,Magic,runID,'.png'];
            SIM21Analysis.plotData(outputPath,figName,XYData,figSettings);
        end
        
        
        function dataMat = getZData(dataPath,outputPath,z,Magic,runID)
            % Create or import mean data matrix
            dataName = [outputPath,Magic,'_Data',runID,'.mat'];
            
            % Check if output exists and load
            if exist(dataName, 'file') == 2
                SIM21Analysis.message('loading data');
                dataMat=importdata(dataName);
            else
                % Calculate Mean
                SIM21Analysis.message('calculating mean');
                dataMat=cat(1,z,zeros(1,length(z)));
                for i = 1:length(dataMat)
                    fileName = SIM21Analysis.genDataFileName(dataPath,Magic,runID,dataMat(1,i));
                    if exist(fileName, 'file') == 2
                        dataMat(2,i) = mean(mean(mean(importdata(fileName))));
                    else % for half run simulations
                        SIM21Analysis.errorMessage(['Error: Missing ',Magic,' z=',num2str(dataMat(1,i)),' file']);
                        dataMat(2,i) = NaN;
                    end
                    if mod(i,10) == 0
                        SIM21Analysis.message(['    ',num2str(i),' / ',num2str(length(dataMat))]);
                    end
                end
                
                % Save Output
                SIM21Analysis.message('saving data');
                save(dataName,'dataMat');
            end
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
        
        
        function plotData(outputPath,outputName,XYData,figSettings)
            % Plot data by XY coordinates
            SIM21Analysis.message('plotting');
            f=figure();
            
            x = XYData(1,:);
            y = XYData(2,:);
            
            % Set logarithmic axes
            if isfield(figSettings,'log')
                if isequal(figSettings.log,'yx')
                    figSettings.log = 'xy';
                end
                switch figSettings.log
                    case 'x'
                        semilogx(x,y);
                    case 'y'
                        semilogy(x,y);
                    case 'xy'
                        loglog(x,y)
                    otherwise
                        plot(x,y);
                end
            else
                plot(x,y);
            end
            hold on;
            
            % Interesting lines (for example zero line)
            if isfield(figSettings,'lines')
                for figLine = figSettings.lines
                    plot(figLine{1}(1,:),figLine{1}(2,:),':k','lineWidth',0.1);
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
            ax = gca;
            if isfield(figSettings,'xTick');
                ax.XTick = figSettings.xTick;
            else
                ax.XTick = floor(min(x)/10)*10:10:max(x);
            end
            if isfield(figSettings,'yTick')
                ax.YTick = figSettings.yTick;
            else
                numOfYTicks = 5;
                yhop = (max(y)-min(y))/numOfYTicks;
                ymag = 10^floor(log10(yhop));
                yhop = floor(yhop/ymag)*ymag;
                ax.YTick = (ceil(min(y) / yhop) * yhop):yhop:max(y);
            end
            
            % Titles and Labels
            if isfield(figSettings,'title');
                title(figSettings.title,'FontSize',18);
            end
            if isfield(figSettings,'xLabel');
                xlabel(figSettings.xLabel,'FontSize',12);
            end
            if isfield(figSettings,'yLabel');
                ylabel(figSettings.yLabel,'FontSize',12);
            end
            
            hold off;
            
            SIM21Analysis.message('saving plot');
            saveas(f,outputName);
        end
        
        
        function specialParams=calcSpecialParams(dataPath,outputPath,runID)
            %Calculate interesting parameters for specific run
            SIM21Analysis.message('calculating parameters');
            interpStep=0.001;
            
            xHIData = SIM21Analysis.interpData(SIM21Analysis.getZData(dataPath,outputPath,SIM21Analysis.xHIZ,SIM21Analysis.xHIMagic,runID),interpStep);
            TKData = SIM21Analysis.interpData(SIM21Analysis.getZData(dataPath,outputPath,SIM21Analysis.TKZ,SIM21Analysis.TKMagic,runID),interpStep);
            T21cmData = SIM21Analysis.interpData(SIM21Analysis.getZData(dataPath,outputPath,SIM21Analysis.T21cmZ,SIM21Analysis.T21cmMagic,runID),interpStep);
            
            % MIN / MAX T21cm
            minT21cmInd = find(T21cmData(2,:)==min(T21cmData(2,:)));
            maxT21cmInd = find(T21cmData(2,:)==max(T21cmData(2,:)));
            specialParams.minT21cm.z = T21cmData(1,minT21cmInd);
            specialParams.minT21cm.T = T21cmData(2,minT21cmInd);
            specialParams.maxT21cm.z = T21cmData(1,maxT21cmInd);
            specialParams.maxT21cm.T = T21cmData(2,maxT21cmInd);
            
            % MIN / MAX Slope
            slope = diff(T21cmData(2,:)) ./ diff(T21cmData(1,:));
            minSlopeInd = find(slope==min(slope));
            maxSlopeInd = find(slope==max(slope));
            specialParams.minSlope.z = T21cmData(1,minSlopeInd);
            specialParams.minSlope.T = T21cmData(2,minSlopeInd);
            specialParams.minSlope.slope = slope(minSlopeInd);
            specialParams.maxSlope.z = T21cmData(1,maxSlopeInd);
            specialParams.maxSlope.T = T21cmData(2,maxSlopeInd);
            specialParams.maxSlope.slope = slope(maxSlopeInd);
            
            % MIN TK
            minTKInd = find(TKData(2,:)==min(TKData(2,:)));
            specialParams.minTK.z = TKData(1,minTKInd);
            specialParams.minTK.T = TKData(2,minTKInd);

            % X Crossing
            xCrossInd = find(diff(sign(T21cmData(2,:))));
            specialParams.xCross.z = T21cmData(1,xCrossInd);
            specialParams.xCross.z1 = max(specialParams.xCross.z); % in case there are several zero crossings
            
            % xHI Percentage
            xHI75Ind = find(diff(sign(xHIData(2,:)-0.75)));
            specialParams.xHI75.z = xHIData(1,xHI75Ind);
            xHI50Ind = find(diff(sign(xHIData(2,:)-0.50)));
            specialParams.xHI50.z = xHIData(1,xHI50Ind);
            xHI25Ind = find(diff(sign(xHIData(2,:)-0.25)));
            specialParams.xHI25.z = xHIData(1,xHI25Ind);
            xHI0Ind = max(find(xHIData(2,:)==0));
            specialParams.xHI0.z = xHIData(1,xHI0Ind);
            
            % Heating Transition: Tcmb = Tk
            TCMBData = SIM21Analysis.interpData(cat(1,SIM21Analysis.TKZ,SIM21Gets.getTcmb(SIM21Analysis.TKZ)),interpStep);
            THTInd = find(diff(sign(TKData(2,:)-TCMBData(2,:))));
            specialParams.THT = TKData(:,THTInd);
            specialParams.THT.z = TKData(1,THTInd);
            specialParams.THT.T = TKData(2,THTInd);
            
            % Save matrix
            fileName = [outputPath,'specialParams',runID];
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
        
        
        function tauCMB = checkTau(dataPath,outputPath,runID)
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

            xHImat = SIM21Analysis.getZData(dataPath,outputPath,SIM21Analysis.xHIZ,SIM21Analysis.xHIMagic,runID);
            xHImat = [zeros(1,60),xHImat(2,:)];
            tauCMB = sum((1-xHImat(1:151)).*ne(1:151).*sigmaT.*dldz(1:151))*0.1+sum((1-xHImat(152:156)).*ne(152:156).*sigmaT.*dldz(152:156))*1;
        end
        
        
        function dataFileName = genDataFileName(dataPath, Magic, runID, z)
            % Get raw data matrix file name
            if runID(1) ~= '_'
                runID = ['_',runID];
            end
            dataFileName = [dataPath,Magic,'_',num2str(z),runID,'.mat'];
        end
        
        
        function message(msg)
            % SIM21Analysis class output
            if false
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
