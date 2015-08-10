classdef analysis
    properties(Constant)
        %Data Matrix settings
        xHIMagic = 'xHI';
        xHIZ=cat(2,6:0.1:14.9,15:60);
        TKMagic = 'TK';
        TKZ=5:64;
        T21cmMagic = 'T21cm';
        T21cmZ=6:60;
        NeutMagic = 'Neut';
        NeutZ=cat(2,6:0.1:14.9,15:60);
    end

    methods(Static)
        function analyse()
            dataPath = '/scratch300/matanlotem/Data/';
            outputPath = '/scratch300/matanlotem/Analysis/';
            runID = '0_0.05_1_16.5_1_1_0.075_0_0_2_1_2';
            
            analysis.plotZGraphs(dataPath,outputPath,runID);
            specialParams=analysis.calcSpecialParams(dataPath,outputPath,runID);
            %analysis.plotPowerSpectrum(dataPath,outputPath,runID);
        end


        function specialParams=calcSpecialParams(dataPath,outputPath,runID)
            %Calculate interesting parameters for specific run
            
            analysis.message('calculating parameters');
            interpStep=0.001;

            xHIData = analysis.interpData(analysis.getZData(dataPath,outputPath,analysis.xHIZ,analysis.xHIMagic,runID),interpStep);
            TKData = analysis.interpData(analysis.getZData(dataPath,outputPath,analysis.TKZ,analysis.TKMagic,runID),interpStep);
            T21cmData = analysis.interpData(analysis.getZData(dataPath,outputPath,analysis.T21cmZ,analysis.T21cmMagic,runID),interpStep);
              
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
            TCMBData = analysis.interpData(cat(1,analysis.TKZ,2.725*(1+analysis.TKZ)),interpStep);
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
            analysis.plotByZ(dataPath,outputPath,analysis.TKZ,analysis.TKMagic,runID,'TK(z)','1+z','TK [mK]'); %TK
            analysis.plotByZ(dataPath,outputPath,analysis.T21cmZ,analysis.T21cmMagic,runID,'T21cm(z)','1+z','T21cm [mK]'); %T21cm
            analysis.plotByZ(dataPath,outputPath,analysis.xHIZ,analysis.xHIMagic,runID,'xHI(z)','1+z','xHi'); %xHI
        end


        function plotByZ(dataPath,outputPath,z,magic,runID,figTitle,figXLabel,figYLabel)
            % Plot different graphs as function of Z
            
            analysis.message(['=== STARTING ',magic,' ===']);
            interpStep = 0.1;
            dataMat = analysis.getZData(dataPath,outputPath,z,magic,runID);
            XYData = analysis.interpData(dataMat,interpStep);
            XYData(1,:) = XYData(1,:) + 1; % z+1
            figName = strcat(outputPath,magic,'_',runID,'.png');
            analysis.plotData(outputPath,figName,XYData,figTitle,figXLabel,figYLabel);
        end


        function dataMat=getZData(dataPath,outputPath,z,magic,runID)
            % Create or import mean data matrix
            
            dataName = strcat(outputPath,magic,'_Data_',runID,'.mat');
            
            % Check if output exists and load
            if exist(dataName, 'file') == 2
                analysis.message('loading data');
                dataMat=importdata(dataName);
            else
                % Calculate Mean
                analysis.message('calculating mean');
                dataMat=cat(1,z,zeros(1,length(z)));
                for i = 1:length(dataMat)
                    dataMat(2,i)=mean(mean(mean(importdata(analysis.genDataFileName(dataPath,magic,runID,dataMat(1,i))))));
                    if mod(i,10) == 0
                        analysis.message(['    ',num2str(i),' / ',num2str(length(dataMat))]);
                    end
                end
                
                % Save Output
                analysis.message('saving data');
                save(dataName,'dataMat');
            end
        end


        function XYData=interpData(dataMat,interpStep)
            % Interpolate data matrix
            analysis.message('interpolating');
            x = dataMat(1,:);
            y = dataMat(2,:);
            x_i = min(x):interpStep:max(x);
            y_i = interp1(x,y,x_i,'spline');
            XYData = cat(1,x_i,y_i);
        end


        function plotData(outputPath,outputName,XYData,figTitle,figXLabel,figYLabel)
            % Plot data by XY coordinates
            analysis.message('plotting');
            f=figure();
            hold on;
            plot(XYData(1,:),XYData(2,:));
            title(figTitle,'FontSize',18);
            xlabel(figXLabel,'FontSize',12);
            ylabel(figYLabel,'FontSize',12);
            hold off;

            analysis.message('saving plot');
            saveas(f,outputName);
        end


        function plotPowerSpectrum(dataPath,outputPath,runID)
            analysis.message('plotting power spectrum');
            
            figName = strcat(outputPath,'PW_',runID,'.png');
            PS = (K(Ind1).^3.*squeeze(real(PowerMat(:,:,Ind1)))/(2*pi^2));
            PSiso = (K(Ind1).^3.*squeeze(real(PowerMat_iso(:,:,Ind1)))/(2*pi^2));
            PSx = (K(Ind1).^3.*squeeze(real(PowerMat_X(:,:,Ind1)))/(2*pi^2) );
            PSdel = (K(Ind1).^3.*squeeze(real(PowerMat_del(:,:,Ind1)))/(2*pi^2));
            zPS = min(zspec):0.01:max(zspec); 

            delx = 0.04;
            x2=min(zspec):delx:max(zspec); 
            Ind = floor(1+delx*(0:length(x2)-1)/0.01);

            f=figure();
            PS1 = interp1(zspec,(PS(2,:)),zPS,'spline');
            hl1 = loglog(1+x2,PS1(Ind).*(x2>6.9),'Color','r','LineWidth',1);
            hold on
            PS1 = interp1(zspec,(PSiso(2,:)),zPS,'spline');
            loglog(1+x2,PS1(Ind).*(x2>6.9),'Color',[0,127/255,0],'LineStyle','-','LineWidth',1);
            PS1 = interp1(zspec,real(PSx(2,:)),zPS,'spline');
            loglog(1+x2,PS1(Ind).*(x2>6.9),'Color','b','LineStyle','-','LineWidth',1);% color green [0,127/255,0]
            PS1 = interp1(zspec,-real(PSx(2,:)),zPS,'spline');
            loglog(1+x2,PS1(Ind).*(x2>6.9),'Color','b','LineStyle',':','LineWidth',2);% color green [0,127/255,0
            PS1 = interp1(zspec,(PSdel(2,:)),zPS,'spline');
            loglog(1+x2,PS1(Ind).*(x2>6.9),'Color','k','LineStyle','-','LineWidth',1);% color green [0,127/255,0]
            xlim([7,40])
            ylim([1e-4,1e3])
            xlabel('1+z','FontSize',12);
            ylabel('k^3P(k)/2\pi^2 [mK^2]','FontSize',12);
            title('Atomic, Old spectrum','FontSize',15,'FontWeight','Bold')
            hold off
            
            
            analysis.message('saving plot');
            saveas(f,outputName);
         end


        function dataFileName = genDataFileName(dataPath, magic, runID, z)
            % Get raw data matrix file name
            dataFileName = strcat(dataPath,magic,'_',num2str(z),'_',runID,'.mat');
        end

        function message(msg)
            % Analysis class output
            disp(msg)
        end

    end
end
