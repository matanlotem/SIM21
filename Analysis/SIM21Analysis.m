classdef SIM21Analysis
    properties(Constant)
        %Data Matrix settings
        xHIMagic = 'xHI';
        xHIZ = [6:0.1:15,16:60];
        TKMagic = 'TK';
        TKZ = 5:64;
        T21cmMagic = 'T21cm';
        T21cmZ = 6:60;
        NeutMagic = 'Neut';
        NeutZ = [6:0.1:15,16:60];
        %PwSpZ = [5:0.1:15,16:40];
        PwSpZ = [6:0.1:15,16:40];
    end
    
    methods(Static)
        function analyse()
            dataPath = '/scratch100/matanlotem/Data/';
            outputPath = '/scratch100/matanlotem/Analysis/';
            runID = '0_0.05_1_16.5_1_1_0.075_0_0_2_1_2';
            
            %SIM21Analysis.plotZGraphs(dataPath,outputPath,runID);
            %specialParams=SIM21Analysis.calcSpecialParams(dataPath,outputPath,runID);
            SIM21Analysis.plotPowerSpectrumByK(outputPath,runID,0.1);
            SIM21Analysis.plotPowerSpectrumByK(outputPath,runID,0.5);
        end
        
        
        function specialParams=calcSpecialParams(dataPath,outputPath,runID)
            %Calculate interesting parameters for specific run
            SIM21Analysis.message('calculating parameters');
            interpStep=0.001;
            
            xHIData = SIM21Analysis.interpData(SIM21Analysis.getZData(dataPath,outputPath,SIM21Analysis.xHIZ,SIM21Analysis.xHIMagic,runID),interpStep);
            TKData = SIM21Analysis.interpData(SIM21Analysis.getZData(dataPath,outputPath,SIM21Analysis.TKZ,SIM21Analysis.TKMagic,runID),interpStep);
            T21cmData = SIM21Analysis.interpData(SIM21Analysis.getZData(dataPath,outputPath,SIM21Analysis.T21cmZ,SIM21Analysis.T21cmMagic,runID),interpStep);
            
            % MIN / MAX
            minT21cmInd = find(T21cmData(2,:)==min(T21cmData(2,:)));
            maxT21cmInd = find(T21cmData(2,:)==max(T21cmData(2,:)));
            specialParams.maxT21cm = T21cmData(:,maxT21cmInd);
            specialParams.minT21cm = T21cmData(:,minT21cmInd);
            
            % MIN / MAX Slope
            slope = diff(T21cmData(2,:)) ./ diff(T21cmData(1,:));
            minSlopeInd = find(slope==min(slope));
            maxSlopeInd = find(slope==max(slope));
            specialParams.minSlope = cat(1,T21cmData(:,minSlopeInd),slope(minSlopeInd));
            specialParams.maxSlope = cat(1,T21cmData(:,maxSlopeInd),slope(maxSlopeInd));
            
            % X Crossing
            xCrossInd = find(diff(sign(T21cmData(2,:))));
            specialParams.xCross = T21cmData(1,xCrossInd);
            
            % xHI Percentage
            xHI75Ind = find(diff(sign(xHIData(2,:)-0.75)));
            specialParams.xHI75 = xHIData(1,xHI75Ind);
            xHI50Ind = find(diff(sign(xHIData(2,:)-0.50)));
            specialParams.xHI50 = xHIData(1,xHI50Ind);
            xHI25Ind = find(diff(sign(xHIData(2,:)-0.25)));
            specialParams.xHI25 = xHIData(1,xHI25Ind);
            
            % Heating Transition: Tcmb = Tk
            TCMBData = SIM21Analysis.interpData(cat(1,SIM21Analysis.TKZ,2.725*(1+SIM21Analysis.TKZ)),interpStep);
            THTInd = find(diff(sign(TKData(2,:)-TCMBData(2,:))));
            specialParams.THT = TKData(:,THTInd);
            
            % Save matrix
            fileName = [outputPath,'specialParams_',runID];
            save([fileName,'.mat'],'specialParams');
            
            % Save readable output
            niceOutput = ['Min T\n\tz = ',num2str(specialParams.minT21cm(1)),'\n\tT = ',num2str(specialParams.minT21cm(2)),'\n',...
                          'Max T\n\tz = ',num2str(specialParams.maxT21cm(1)),'\n\tT = ',num2str(specialParams.maxT21cm(2)),'\n',... 
                          'Min Slope\n\tz = ',num2str(specialParams.minSlope(1)),'\n\tT = ',num2str(specialParams.minSlope(2)),'\n\tSlope = ',num2str(specialParams.minSlope(3)),'\n',... 
                          'Max Slope\n\tz = ',num2str(specialParams.maxSlope(1)),'\n\tT = ',num2str(specialParams.maxSlope(2)),'\n\tSlope = ',num2str(specialParams.maxSlope(3)),'\n',...
                          '0 Crossing\n\tz = ',num2str(specialParams.xCross),'\n',...
                          'xHI Percentage\n\t75%% z = ',num2str(specialParams.xHI75),'\n\t50%% z = ',num2str(specialParams.xHI50),'\n\t25%% z = ',num2str(specialParams.xHI25),'\n',...
                          'Heating Transition\n\tz = ',num2str(specialParams.THT(1)),'\n\tT = ',num2str(specialParams.THT(2)),'\n'];
            
            fid = fopen([fileName,'.txt'],'w');
            fwrite(fid,sprintf(niceOutput));
            fclose(fid);
        end
        
        
        function plotZGraphs(dataPath,outputPath,runID)
            % Plot TK, T21cm and xHI Graphs
            SIM21Analysis.plotByZ(dataPath,outputPath,SIM21Analysis.TKZ,SIM21Analysis.TKMagic,runID,'TK(z)','1+z','TK [mK]'); %TK
            SIM21Analysis.plotByZ(dataPath,outputPath,SIM21Analysis.T21cmZ,SIM21Analysis.T21cmMagic,runID,'T21cm(z)','1+z','T21cm [mK]'); %T21cm
            SIM21Analysis.plotByZ(dataPath,outputPath,SIM21Analysis.xHIZ,SIM21Analysis.xHIMagic,runID,'xHI(z)','1+z','xHi'); %xHI
        end
        
        
        function plotByZ(dataPath,outputPath,z,magic,runID,figTitle,figXLabel,figYLabel)
            % Plot different graphs as function of Z
            
            SIM21Analysis.message(['=== STARTING ',magic,' ===']);
            interpStep = 0.1;
            dataMat = SIM21Analysis.getZData(dataPath,outputPath,z,magic,runID);
            XYData = SIM21Analysis.interpData(dataMat,interpStep);
            XYData(1,:) = XYData(1,:) + 1; % z+1
            figName = strcat(outputPath,magic,'_',runID,'.png');
            SIM21Analysis.plotData(outputPath,figName,XYData,figTitle,figXLabel,figYLabel);
        end
        
        
        function dataMat=getZData(dataPath,outputPath,z,magic,runID)
            % Create or import mean data matrix
            
            dataName = strcat(outputPath,magic,'_Data_',runID,'.mat');
            
            % Check if output exists and load
            if exist(dataName, 'file') == 2
                SIM21Analysis.message('loading data');
                dataMat=importdata(dataName);
            else
                % Calculate Mean
                SIM21Analysis.message('calculating mean');
                dataMat=cat(1,z,zeros(1,length(z)));
                for i = 1:length(dataMat)
                    dataMat(2,i)=mean(mean(mean(importdata(SIM21Analysis.genDataFileName(dataPath,magic,runID,dataMat(1,i))))));
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
        
        
        function plotData(outputPath,outputName,XYData,figTitle,figXLabel,figYLabel)
            % Plot data by XY coordinates
            SIM21Analysis.message('plotting');
            f=figure();
            hold on;
            plot(XYData(1,:),XYData(2,:));
            title(figTitle,'FontSize',18);
            xlabel(figXLabel,'FontSize',12);
            ylabel(figYLabel,'FontSize',12);
            hold off;
            
            SIM21Analysis.message('saving plot');
            saveas(f,outputName);
        end
        
        
        function plotPowerSpectrumByK(outputPath,runID,k)
            SIM21Analysis.message('plotting power spectrum');
            
            figName = [outputPath,'PwSpK_',num2str(k),'_',runID,'.png'];
            MK = SIM21Analysis.importMatrix('K');
            KInd = find(diff(sign(MK-k)));
            
            delx = 0.04;
            x2=min(SIM21Analysis.PwSpZ):delx:max(SIM21Analysis.PwSpZ); 
            Ind = floor(1+delx*(0:length(x2)-1)/0.01);
            zPS = min(SIM21Analysis.PwSpZ):0.01:max(SIM21Analysis.PwSpZ); 
            
            f=figure();
            hold on
            
            function doPowerSpectrum(MName,Msign,lineColor,lineStyle,lineWidth)
                PowerMat = importdata([outputPath,MName,'_',runID,'.mat']);
                PS = (MK(KInd).^3.*squeeze(real(PowerMat(:,KInd)))/(2*pi^2));
                PS1 = interp1(SIM21Analysis.PwSpZ,Msign*PS',zPS,'spline');
                loglog(1+x2,PS1(Ind).*(x2>6.9),'Color',lineColor,'LineStyle',lineStyle,'LineWidth',lineWidth);
            end
            
            doPowerSpectrum('PowerMat',1,'r','-',1);
            doPowerSpectrum('PowerMat_iso',1,[0,127/255,0],'-',1);
            doPowerSpectrum('PowerMat_X',1,'b','-',1);
            doPowerSpectrum('PowerMat_X',-1,'b',':',2);
            doPowerSpectrum('PowerMat_del',1,'k','-',1);
            
%            PowerMat = importdata([outputPath,'PowerMat_',runID,'.mat']);
%            PowerMat_iso = importdata([outputPath,'PowerMat_iso_',runID,'.mat']);
%            PowerMat_X = importdata([outputPath,'PowerMat_X_',runID,'.mat']);
%            PowerMat_del = importdata([outputPath,'PowerMat_del_',runID,'.mat']);
%            
%            PS = (MK(KInd).^3.*squeeze(real(PowerMat(:,KInd)))/(2*pi^2));
%            PSiso = (MK(KInd).^3.*squeeze(real(PowerMat_iso(:,KInd)))/(2*pi^2));
%            PSx = (MK(KInd).^3.*squeeze(real(PowerMat_X(:,KInd)))/(2*pi^2) );
%            PSdel = (MK(KInd).^3.*squeeze(real(PowerMat_del(:,KInd)))/(2*pi^2));
%            
%            PS1 = interp1(SIM21Analysis.PwSpZ,(PS'),zPS,'spline');
%            hl1 = loglog(1+x2,PS1(Ind).*(x2>6.9),'Color','r','LineWidth',1);
%            
%            PS1 = interp1(SIM21Analysis.PwSpZ,(PSiso'),zPS,'spline');
%            loglog(1+x2,PS1(Ind).*(x2>6.9),'Color',[0,127/255,0],'LineStyle','-','LineWidth',1);
%            
%            PS1 = interp1(SIM21Analysis.PwSpZ,real(PSx'),zPS,'spline');
%            loglog(1+x2,PS1(Ind).*(x2>6.9),'Color','b','LineStyle','-','LineWidth',1);% color green [0,127/255,0]
%            PS1 = interp1(SIM21Analysis.PwSpZ,-real(PSx'),zPS,'spline');
%            loglog(1+x2,PS1(Ind).*(x2>6.9),'Color','b','LineStyle',':','LineWidth',2);% color green [0,127/255,0
%            PS1 = interp1(SIM21Analysis.PwSpZ,(PSdel'),zPS,'spline');
%            loglog(1+x2,PS1(Ind).*(x2>6.9),'Color','k','LineStyle','-','LineWidth',1);% color green [0,127/255,0]
            xlim([7,40])
            ylim([1e-4,1e3])
            xlabel('1+z','FontSize',12);
            ylabel('k^3P(k)/2\pi^2 [mK^2]','FontSize',12);
            title('Atomic, Old spectrum','FontSize',15,'FontWeight','Bold')
            hold off
            
            SIM21Analysis.message('saving plot');
            saveas(f,figName);
         end
        
        
        function dataFileName = genDataFileName(dataPath, magic, runID, z)
            % Get raw data matrix file name
            dataFileName = strcat(dataPath,magic,'_',num2str(z),'_',runID,'.mat');
        end
        
        function M = importMatrix(MName)
            M = importdata(['Matrices/',MName,'.mat']);
        end
        
        function message(msg)
            % Analysis class output
            disp(msg)
        end
    end
end
