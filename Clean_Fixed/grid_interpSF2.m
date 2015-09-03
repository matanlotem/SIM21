function [fgas,fgasp,fgasMQ,fgasMQp] = grid_interpSF2(flag,flagM,JLW21,z,Ispec,fstarM,fstarA,FSfunc,photoheatingOn,photoheatingVersion,calcfgasp)
    global pathname_Data2
    global pathname_DataBackgrounds
    
    V1zINchanges=1;

    %%% MATAN CHANGE - 2015/09/02
    %MQ_Vc=16.5; %MQ doesn't form below atomic cooling mass according to Tal Alexander's work.
    MQ_Vc=max(16.5,flagM);
    %%% END CHANGE
    global grid_interpSF2_returnOnlyFgas; % <-- REMOVE DEFINITION BELOW in previous versions

    if isempty(grid_interpSF2_returnOnlyFgas)
        grid_interpSF2_returnOnlyFgas = 0;
    end

    % if 1
    % persistent calculate_SF2;
    % if isempty(calculate_SF2)
    %     calculate_SF2 = 0;
    % end
    % 
    % persistent data_list;
    % if isempty(data_list)
    %     data_list = containers.Map;
    % end
    % 
    % persistent zer;
    % if isempty(zer)
    %     N = size(JLW21,1);
    %     zer = zeros(N,N,N);
    % else
    %     if ~isequal(size(zer),size(JLW21))
    %         N = size(JLW21,1);
    %         zer = zeros(N,N,N);
    %     end
    % end
    % 
    % if (calculate_SF2==0) && (isequal(zer,JLW21))
    %     key = [num2str(grid_interpSF2_returnOnlyFgas) '_' num2str(size(JLW21)) '_' num2str(flag) '_' num2str(flagM) '_' num2str(z) '_' num2str(fstarM) '_' num2str(fstarA)];
    %     if data_list.isKey(key)
    %         t = data_list(key);
    %         fgas = t(:,:,:,1);
    %         fgasp = t(:,:,:,2);
    %     else
    %         calculate_SF2 = 1;
    %         [fgas,fgasp,fgasMQ,fgasMQp] = grid_interpSF2(flag,flagM,JLW21,z,Ispec,fstarM,fstarA,FSfunc,photoheatingOn,photoheatingVersion,1);
    %         N = size(JLW21,1);
    %         v = zeros(N,N,N,2);
    %         v(:,:,:,1) = fgas;
    %         if ~isequal(fgasp,[])
    %             v(:,:,:,2) = fgasp;
    %         end
    %         data_list(key) = v;
    %     end
    %     grid_interpSF2_returnOnlyFgas = 0;
    %     calculate_SF2 = 0;
    %     return;
    % end
    % end


    global vbc_cube
    global delta_cube

    if (z>15)
        delmax = 1.0;
        Nres = 100;
        deltas = linspace(-delmax,delmax,Nres);
        vbc = linspace(0,3.8,Nres);
        NresM = 45;
        if (z<66)
            Mass = logspace(log10(1e5),log10(1e11),NresM);%Mass = logspace(log10(1e5),log10(7e8),NresM);
        else
            Mass = logspace(log10(1e5),log10(1e8),NresM);
        end
    else
         
        delmax = 1.0;
        Nres = 100;
        deltas = linspace(-delmax,2*delmax,3*Nres/2);
        vbc = linspace(0,3.8,Nres);
        NresM = 45;
        Mass = logspace(log10(1e5),log10(1e11),NresM);% Mass = logspace(log10(1e5),log10(1e10),NresM);
    end

    delta_c=1.686-1.686*(1+z)/2000;
    dz = 0.0001;
    g = LWgetDz(z)/LWgetDz(40);%LWgetDz(40) = 0.031880213
    gp = LWgetDz(z+dz)/LWgetDz(40);

    % --- change 09.19.2013 -----------------------%
    %McubeGM = LWgetMcubeGM(flag ,z, JLW21,flagM); 
    McubeGM = LWgetMcubeVc(flag, z, JLW21,flagM);% here flagM is the value of circular velocity
    McubeGM_MQ = LWgetMcubeVc(flag, z, JLW21,MQ_Vc);
    %----------------------------------------------%


    global DEBUG;
    if isempty(DEBUG)
        DEBUG = 0;
    end


    %persistent directory;
    %if isempty(directory)
    %    directory = '/scratch300/matanlotem/Data/';
    %end
    directory = pathname_Data2;
    persistent filePrefix;
    if isempty(filePrefix)
        filePrefix = 'xHI_';
    end
    global gridInterpKey;
    if isempty(gridInterpKey)
        photoheatingOn = 0;
    end
    photoheating=0;
    if photoheatingOn==1
        photoheating = photoheatingVersion;
    end

    zMax = 60; % Not sure if 60 is the right place to start
    zMin = 6;

    N = size(McubeGM,1);

    persistent McritHistory;
    if isempty(McritHistory)
        McritHistory = containers.Map(double(1), zeros(size(McubeGM))); remove(McritHistory,1);
    end

    persistent calculateForAllZ;
    if isempty(calculateForAllZ)
        if z==zMax
            calculateForAllZ = 0;
        else
            calculateForAllZ = 1;
        end
    end
    calculateForAllZ = 0;
    zMaxFileExists = exist([directory filePrefix num2str(zMax) gridInterpKey],'file')==2;

    if photoheating==1 && zMaxFileExists


        skipMcrit = McritHistory.isKey(z);
        skipMcrit=0;
        if skipMcrit
            Mcrit = McritHistory(z);
        end

        if ~skipMcrit && z<zMax
            

            if V1zINchanges==0
                Mcrit = calculateMcritV1fixed(directory,filePrefix, zMax,N,z);
            elseif V1zINchanges==1
                Mcrit = calculateMcritV1changes(directory,filePrefix, zMax,N,z);
            end
            
            McritHistory(z) = Mcrit;
        end

        if z<zMax
            %% Set min mass as max(Mcool,Mcrit)
            McubeGM = max(McubeGM,Mcrit);
            McubeGM_MQ = max(McubeGM_MQ,Mcrit);
        end
    end

    if(z>65)%||z==15)
        gridA=0;gridM=0;
        gridAMQ=0;gridMMQ=0;
        load(strcat(pathname_DataBackgrounds,'gridM_',num2str(z),'.mat'));
        load(strcat(pathname_DataBackgrounds,'gridA_',num2str(z),'.mat'));
        if Ispec>2
            load(strcat(pathname_DataBackgrounds,'gridMMQ_',num2str(z),'.mat')); %Q
            load(strcat(pathname_DataBackgrounds,'gridAMQ_',num2str(z),'.mat'));   %Q
            %gridMMQ=zeros(size(gridAMQ)); %MQ doesn't form below atomic cooling mass according to Tal Alexander's work
        end
    else
      
        if FSfunc==1
            gridA=0;gridM=0;    
            gridAMQ=0;gridMMQ=0;
            load(strcat(pathname_DataBackgrounds,'gridM10_',num2str(z),'.mat'));
            load(strcat(pathname_DataBackgrounds,'gridA10_',num2str(z),'.mat'));
            if Ispec>2
                load(strcat(pathname_DataBackgrounds,'gridMMQ10_',num2str(z),'.mat')); %Q
                load(strcat(pathname_DataBackgrounds,'gridAMQ10_',num2str(z),'.mat'));   %Q
                %gridMMQ=zeros(size(gridAMQ)); %MQ doesn't form below atomic cooling mass according to Tal Alexander's work
            end
        end
        if FSfunc==2
            gridA=0;gridM=0;
            gridAMQ=0;gridMMQ=0;
            load(strcat(pathname_DataBackgrounds,'gridM10Sharp_',num2str(z),'.mat'));
            load(strcat(pathname_DataBackgrounds,'gridA10Sharp_',num2str(z),'.mat')); 
            if Ispec>2
                load(strcat(pathname_DataBackgrounds,'gridM10SharpMQ_',num2str(z),'.mat')); %Q
                load(strcat(pathname_DataBackgrounds,'gridA10SharpMQ_',num2str(z),'.mat')); %Q
                %gridMMQ=zeros(size(gridAMQ)); %MQ doesn't form below atomic cooling mass according to Tal Alexander's work
            end  
        end
    end


     
    fgasA = exp(interp3(vbc,deltas,Mass,log(gridA),(flag+1e-10)*vbc_cube,...
        max(min(delta_cube*g,0.99*delta_c*ones(size(delta_cube))),-0.99*ones(size(delta_cube))),McubeGM,'linear'));

    fgasM = exp(interp3(vbc,deltas,Mass,log(gridM),(flag+1e-10)*vbc_cube,...
        max(min(delta_cube*g,0.99*delta_c*ones(size(delta_cube))),-0.99*ones(size(delta_cube))),McubeGM,'linear'));
     
    fgas = fstarM*fgasM+fstarA*fgasA;

    fgasMQ=[];
    if Ispec>2
        fgasAMQ = exp(interp3(vbc,deltas,Mass,log(gridAMQ),(flag+1e-10)*vbc_cube,...
            max(min(delta_cube*g,0.99*delta_c*ones(size(delta_cube))),-0.99*ones(size(delta_cube))),McubeGM_MQ,'linear'));%Q

        fgasMMQ = exp(interp3(vbc,deltas,Mass,log(gridMMQ),(flag+1e-10)*vbc_cube,...
            max(min(delta_cube*g,0.99*delta_c*ones(size(delta_cube))),-0.99*ones(size(delta_cube))),McubeGM_MQ,'linear'));%Q

        fgasMQ = fstarM*fgasMMQ+fstarA*fgasAMQ;%Q
    end


    if photoheating==2
       

        zMaxFileExists = exist([directory filePrefix num2str(zMax) gridInterpKey],'file')==2;

        if z<(zMax-1) && zMaxFileExists

            skipMcrit = McritHistory.isKey(z);
            skipMcrit = 0;
            if skipMcrit
                Mcrit = McritHistory(z);
            
                %% Load current ionization
                xHI = loadMatOnce([directory filePrefix num2str(z+1) gridInterpKey],'xHI');
                currentIonization = 1-xHI;
            else
                
                %%
                %% Calculate Mcrit
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [Mcrit,currentIonization] = calculateMcritV2(flag,flagM,JLW21,fstarM,fstarA,FSfunc,directory,filePrefix,zMax,zMin,N,z);

                %%
                %% Save Mcrit
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                McritHistory(z) = Mcrit;
            end
        
            %%
            %% Calculate new Mmin
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            

            
            oldMcubeGM = McubeGM;
            McubeGM = max(McubeGM,Mcrit);
            newMcubeGM = McubeGM;
            
            oldMcubeGM_MQ = McubeGM_MQ;
            McubeGM_MQ = max(McubeGM_MQ,Mcrit);
            newMcubeGM_MQ = McubeGM_MQ;
            
            %% Calculate new fgas
            fgasA = exp(interp3(vbc,deltas,Mass,log(gridA),(flag+1e-10)*vbc_cube,...
                max(min(delta_cube*g,0.99*delta_c*ones(size(delta_cube))),-0.99*ones(size(delta_cube))),McubeGM,'linear'));
            gridA=0;
            fgasM = exp(interp3(vbc,deltas,Mass,log(gridM),(flag+1e-10)*vbc_cube,...
                max(min(delta_cube*g,0.99*delta_c*ones(size(delta_cube))),-0.99*ones(size(delta_cube))),McubeGM,'linear'));
            gridM=0;
            fgasNew = fstarM*fgasM+fstarA*fgasA;

            fgas = fgas.*(1-currentIonization) + fgasNew.*currentIonization;
            
            if Ispec>2

                    
                fgasAMQ = exp(interp3(vbc,deltas,Mass,log(gridAMQ),(flag+1e-10)*vbc_cube,...
                   max(min(delta_cube*g,0.99*delta_c*ones(size(delta_cube))),-0.99*ones(size(delta_cube))),McubeGM_MQ,'linear'));%Q
                gridAMQ=0;%Q
                fgasMMQ = exp(interp3(vbc,deltas,Mass,log(gridMMQ),(flag+1e-10)*vbc_cube,...
                   max(min(delta_cube*g,0.99*delta_c*ones(size(delta_cube))),-0.99*ones(size(delta_cube))),McubeGM_MQ,'linear'));%Q
                gridMMQ=0; %Q
                fgasMQNew = fstarM*fgasMMQ+fstarA*fgasAMQ;%Q
                fgasMQ = fgasMQ.*(1-currentIonization) + fgasMQNew.*currentIonization;
            end
        end
    end


    fgasp=[];
    if calcfgasp==1
        if(z>65)%||z==15)
            gridA=0;gridM=0;
            gridAMQ=0;gridMMQ=0;
            load(strcat(pathname_DataBackgrounds,'gridM_',num2str(z),'_1.mat'));
            load(strcat(pathname_DataBackgrounds,'gridA_',num2str(z),'_1.mat'));
            if Ispec>2
                load(strcat(pathname_DataBackgrounds,'gridMMQ_',num2str(z),'_1.mat')); %Q
                load(strcat(pathname_DataBackgrounds,'gridAMQ_',num2str(z),'_1.mat'));   %Q
                %gridMMQ=zeros(size(gridAMQ)); %MQ doesn't form below atomic cooling mass according to Tal Alexander's work
            end
        else
            if FSfunc==1
                gridA=0;gridM=0;
                gridAMQ=0;gridMMQ=0;
                load(strcat(pathname_DataBackgrounds,'gridM10_',num2str(z),'_1.mat'));
                load(strcat(pathname_DataBackgrounds,'gridA10_',num2str(z),'_1.mat'));
                if Ispec>2
                    load(strcat(pathname_DataBackgrounds,'gridMMQ10_',num2str(z),'_1.mat')); %Q
                    load(strcat(pathname_DataBackgrounds,'gridAMQ10_',num2str(z),'_1.mat'));   %Q
                    %gridMMQ=zeros(size(gridAMQ)); %MQ doesn't form below atomic cooling mass according to Tal Alexander's work
                end
            end
            if FSfunc==2
                gridA=0;gridM=0;
                gridAMQ=0;gridMMQ=0;
                load(strcat(pathname_DataBackgrounds,'gridM10Sharp_',num2str(z),'_1.mat'));
                load(strcat(pathname_DataBackgrounds,'gridA10Sharp_',num2str(z),'_1.mat'));
                
                if Ispec>2
                    load(strcat(pathname_DataBackgrounds,'gridM10SharpMQ_',num2str(z),'_1.mat')); %Q
                    load(strcat(pathname_DataBackgrounds,'gridA10SharpMQ_',num2str(z),'_1.mat')); %Q
                    %gridMMQ=zeros(size(gridAMQ)); %MQ doesn't form below atomic cooling mass according to Tal Alexander's work
                end
            end
        end


        if photoheating==2
        
            if z<(zMax-1) && zMaxFileExists
                McubeGM = oldMcubeGM;
                McubeGM_MQ = oldMcubeGM_MQ;
            end
        end

        fgasAp = exp(interp3(vbc,deltas,Mass,log(gridA),(flag+1e-10)*vbc_cube,...
            max(min(delta_cube*gp,0.99*delta_c*ones(size(delta_cube))),-0.99*ones(size(delta_cube))),McubeGM,'linear'));
        fgasMp = exp(interp3(vbc,deltas,Mass,log(gridM),(flag+1e-10)*vbc_cube,...
            max(min(delta_cube*gp,0.99*delta_c*ones(size(delta_cube))),-0.99*ones(size(delta_cube))),McubeGM,'linear'));
        fgasp = fstarM*fgasMp+fstarA*fgasAp;
        grid = [];

        fgasMQp=[];
        if Ispec>2
            fgasAMQp = exp(interp3(vbc,deltas,Mass,log(gridAMQ),(flag+1e-10)*vbc_cube,...
                max(min(delta_cube*gp,0.99*delta_c*ones(size(delta_cube))),-0.99*ones(size(delta_cube))),McubeGM_MQ,'linear'));%Q
            fgasMMQp = exp(interp3(vbc,deltas,Mass,log(gridMMQ),(flag+1e-10)*vbc_cube,...
                max(min(delta_cube*gp,0.99*delta_c*ones(size(delta_cube))),-0.99*ones(size(delta_cube))),McubeGM_MQ,'linear'));%Q
            fgasMQp = fstarM*fgasMMQp+fstarA*fgasAMQp; %Q
            McubeGM  = [];
        end
     
        if photoheating==2
            if z<(zMax-1) && zMaxFileExists

                McubeGM = newMcubeGM;
                McubeGM_MQ = newMcubeGM_MQ;

                %% Calculate new fgas
                fgasAp = exp(interp3(vbc,deltas,Mass,log(gridA),(flag+1e-10)*vbc_cube,...
                max(min(delta_cube*gp,0.99*delta_c*ones(size(delta_cube))),-0.99*ones(size(delta_cube))),McubeGM,'linear'));
                fgasMp = exp(interp3(vbc,deltas,Mass,log(gridM),(flag+1e-10)*vbc_cube,...
                max(min(delta_cube*gp,0.99*delta_c*ones(size(delta_cube))),-0.99*ones(size(delta_cube))),McubeGM,'linear'));
                fgaspNew = fstarM*fgasMp+fstarA*fgasAp;

                %%
                %% Update fgasp using the new fgasp with Mcrit
                %% Notes:
                %% + Make sure we keep the same xHI that defined currentIonization
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                fgasp = fgasp.*(1-currentIonization) + fgaspNew.*currentIonization;
                
                
                if Ispec>2
                    fgasAMQp = exp(interp3(vbc,deltas,Mass,log(gridAMQ),(flag+1e-10)*vbc_cube,...
                       max(min(delta_cube*g,0.99*delta_c*ones(size(delta_cube))),-0.99*ones(size(delta_cube))),McubeGM_MQ,'linear'));%Q
                    gridAMQp=0;%Q
                    fgasMMQp = exp(interp3(vbc,deltas,Mass,log(gridMMQ),(flag+1e-10)*vbc_cube,...
                       max(min(delta_cube*g,0.99*delta_c*ones(size(delta_cube))),-0.99*ones(size(delta_cube))),McubeGM_MQ,'linear'));%Q
                    gridMMQp=0; %Q
                    fgasMQpNew = fstarM*fgasMMQp+fstarA*fgasAMQp;%Q
                fgasMQp = fgasMQp.*(1-currentIonization) + fgasMQpNew.*currentIonization;
                end
            end
        end
    end
end
