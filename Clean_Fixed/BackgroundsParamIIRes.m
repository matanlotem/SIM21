
function xHImean = BackgroundsParamIIRes(zcenter,ncube,fstar,flag,flagM,XeffTerm,Ispec,Reion,zeta,feedback,p,pop,FSfunc,photoheatingOn,photoheatingVersion)
    global pathname_Data2
    global pathname_Data1
    load(SIM21Utils.matrixPath('Planck_parameters'));
    %Lpix=3;
    zMAX2 = 50;% 60
    Threshold = 1/zeta;
    global delta_cube
    N=length(delta_cube);

    %zcenter=9.3
    % flagM=3
    % XeffTerm=1
    % Ispec=1
    % Reion=1
    % ncube=0



    %N = 128;
    %JLW21 = zeros(N,N,N);
    % if (feedback==1)
    %         zii =  zcenter;
    %         z0 = (1+zii)/(min(max(p,0.05),1))^(2/3)-1;
    %         J_interp = zeros(2,N,N,N);
    %         zfeed = (ceil(z0)-1:ceil(z0));
    %         % load the LW to calculate the feedback
    %         for jj=1:2    
    %             load(strcat('/scratch/prun/P_DATA/JLW_',num2str(zfeed(jj)),'_',num2str(ncube),'_',num2str(fstar),'_',num2str(flag),...
    %                                 '_',num2str(flagM),'_',num2str(XeffTerm),'_',num2str(Ispec),'_',num2str(Reion),'_',num2str(feedback),'_',num2str(p),'.mat'));% total LyA flux
    %             J_interp(jj,:,:,:) = JLW21; 
    %             JLW21 = [];
    %         end      
    %         JLW21 = squeeze(exp(interp1(zfeed,log(abs(J_interp)),z0))); 
    % else
    %         JLW21 = zeros(N,N,N);
    % end






    Rset = [linspace(1e-10,500/3,110) logspace(log10(550/3),log10(15000/3),20)];
    Rmin = Rset(1:end-1);
    Rmax = Rset(2:end);
    Nshells = length(Rmin);
    if(length(Rmax) ~= length(Rmin))
        pause;
    end

    %zfgas = 5:100;
    zgas = (6:1:16);
     

    Ind1 = find(zgas>zcenter,1,'first');  
    zint =  zgas(Ind1-1:Ind1);  

    %----
    if (length(zint)==2)
        F_interp = zeros(2,N,N,N);
        for iz=1:2
            
            %---
            
            if (feedback==1)
                zii =  zint(iz);
                z0 = (1+zii)/(min(max(p,0.05),1))^(2/3)-1;
                J_interp = zeros(2,N,N,N);
                zfeed = (ceil(z0)-1:ceil(z0));

                % load the LW to calculate the feedback
                for jj=1:2    

                    load(strcat(pathname_Data1,'JLW_',num2str(zfeed(jj)),'_',num2str(ncube),'_',num2str(fstar),'_',num2str(flag),...
                        '_',num2str(flagM),'_',num2str(XeffTerm),'_',num2str(Ispec),'_',num2str(Reion),'_',num2str(feedback),...
                        '_',num2str(p),'_',num2str(pop),'_',num2str(FSfunc),'_',num2str(photoheatingVersion),'.mat'));% total LyA flux
                    J_interp(jj,:,:,:) = JLW21; 
                    JLW21 = [];
                end      
                JLW21 = squeeze(exp(interp1(zfeed,log(abs(J_interp)),z0))); 
                
            else
                JLW21 = zeros(N,N,N);
            end
                        
                        %---
                     
            [fgas_z,fgas_zp] = grid_interpSF2(flag,flagM,feedback*JLW21,zint(iz),Ispec,fstar,fstar,FSfunc,photoheatingOn,photoheatingVersion,0);% fstar = 1 gives collapsed fraction
            F_interp(iz,:,:,:) = fgas_z/fstar;
            fgas_zp=[];        
        end
        fgas = squeeze(exp(interp1(log(1+zint),log(abs(F_interp)+1e-16),log(1+zcenter)))); 
        
         F_interp=[];  
    else
        %---
        if(feedback==1)
            zii =  zint;
            z0 = (1+zii)/(min(max(p,0.05),1))^(2/3)-1;
            J_interp = zeros(2,N,N,N);
            zfeed = (ceil(z0)-1:ceil(z0));

            % load the LW to calculate the feedback
            for jj=1:2    

                load(strcat(pathname_Data1,'JLW_',num2str(zfeed(jj)),'_',num2str(ncube),'_',num2str(fstar),'_',num2str(flag),...
                    '_',num2str(flagM),'_',num2str(XeffTerm),'_',num2str(Ispec),'_',num2str(Reion),'_',num2str(feedback),...
                    '_',num2str(p),'_',num2str(pop),'_',num2str(FSfunc),'_',num2str(photoheatingVersion),'.mat'));% total LyA flux
                J_interp(jj,:,:,:) = JLW21; 
                JLW21 = [];
            end      
            JLW21 = squeeze(exp(interp1(zfeed,log(abs(J_interp)),z0))); 
        else
           JLW21 = zeros(N,N,N);
        end
        %----
        [fgas,fgas_zp] = grid_interpSF2(flag,flagM,feedback*JLW21,zion(indzz),Ispec,fstar,fstar,FSfunc,photoheatingOn,photoheatingVersion,0);% fstar = 1 gives collapsed fraction
        fgas=fgas/fstar;
        fgas_zp=[];
    end
    %----


    % if (length(zint)==2)
    %     F_interp = zeros(2,N,N,N);
    %     for iz=1:2
    %         [fgas_z,fgas_zp] = grid_interpSF2(flag,flagM,feedback*JLW21,zint(iz),1,1);% fstar inside  
    %         F_interp(iz,:,:,:) = fgas_z;
    %         fgas_zp=[];        
    %     end
    %     fgas = squeeze(exp(interp1(log(1+zint),log(abs(F_interp)),log(1+zcenter)))); 
    %     F_interp=[];  
    % else
    %     [fgas,fgas_zp] = grid_interpSF2(flag,flagM,feedback*JLW21,zcenter,1,1);% fstar inside 
    %     fgas_zp=[];
    % end
           
    fcoll_matrix =fftn(fgas);
    Maxfcoll = zeros(N,N,N);
    for ii= 1:Nshells  %integral  
        R=Lpix*(Rmin(ii)+Rmax(ii))/2; % comoving radius of ring
        if(R<70) 
            fcoll = LWgetGasShell2(fcoll_matrix,0,Rmax(ii))/(4*pi*Rmax(ii)^3/3);
            Maxfcoll = max(fcoll,Maxfcoll); 
        end
    end
    fcoll_matrix=[];

    Maxfcoll = max(fgas,Maxfcoll);% adding the pixel
    Neut = (Maxfcoll<Threshold);% 1 if neutral, 0 if fully ionized
    xHI = max(0, (1-zeta*fgas).*Neut);% ionized fraction of each pixel
    xHImean = mean(mean(mean(xHI)));
    save(strcat(pathname_Data2,'xHI_',num2str(zcenter),'_',num2str(ncube),'_',num2str(fstar),'_',num2str(flag),...
                            '_',num2str(flagM),'_',num2str(XeffTerm),'_',num2str(Ispec),'_',num2str(Reion),'_',num2str(feedback),...
                            '_',num2str(p),'_',num2str(pop),'_',num2str(FSfunc),'_',num2str(photoheatingVersion),'.mat'),'xHI');

    save(strcat(pathname_Data2,'Neut_',num2str(zcenter),'_',num2str(ncube),'_',num2str(fstar),'_',num2str(flag),...
                        '_',num2str(flagM),'_',num2str(XeffTerm),'_',num2str(Ispec),'_',num2str(Reion),'_',num2str(feedback),...
                        '_',num2str(p),'_',num2str(pop),'_',num2str(FSfunc),'_',num2str(photoheatingVersion),'.mat'),'Neut');
        
end