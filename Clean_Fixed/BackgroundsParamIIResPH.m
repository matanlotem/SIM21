function xHImean = BackgroundsParamIIResPH(zcenter,ncube,fstar,flag,flagM,XeffTerm,Ispec,Reion,zeta,feedback,p,pop,FSfunc,photoheatingOn,photoheatingVersion)
    global pathname_Data1
    global pathname_Data2
    global ID
    global delta_cube
    N=length(delta_cube);

    load(SIM21Utils.getMatrixPath('Planck_parameters'));
    %%% MATAN CHANGE - 2015/10/06
    xe_grid =[10^(-4),10^(-3.3), 10^(-3),10^(-2.6),10^(-2.3),10^(-2),10^(-1.6),10^(-1.3), 10^(-1),0.5,0.9,0.99,1];% mean of the central pixel
    %%% END CHANGE
    zMAX2 = 60;
    Threshold = 1/zeta;

    Rset = [linspace(1e-10,500/3,110) logspace(log10(550/3),log10(15000/3),20)];
    Rmin = Rset(1:end-1);
    Rmax = Rset(2:end);
    Nshells = length(Rmin);
    if(length(Rmax) ~= length(Rmin))
        pause;
    end

    zgas = (6:1:16);

    if(feedback)
        z0 = (1+zcenter)/(min(max(p,0.05),1))^(2/3)-1;
        J_interp = zeros(2,N,N,N);
        zint = (ceil(z0)-1:ceil(z0));

        % load the LW to calculate the feedback
        for jj=1:2    
            J_interp(jj,:,:,:) = importdata([pathname_Data1,'JLW_',num2str(zint(jj)),ID,'.mat']); 
        end      
        JLW21 = squeeze(exp(interp1(zint,log(abs(J_interp)),z0))); 
    else
        JLW21 = zeros(N,N,N);
    end
    [fgas,fgas_zp] = grid_interpSF2(flag,flagM,feedback*JLW21,zcenter,Ispec,fstar,fstar,FSfunc,photoheatingOn,photoheatingVersion,zeta,0);% fstar = 1 gives collapsed fraction

    fgas=fgas/fstar;
    fgas_zp=[];
   
    fcoll_matrix = fftn(fgas);
    Maxfcoll = zeros(N,N,N);
    %%% MATAN CHANGE - 2015/10/06
    zint = [floor(zcenter),floor(zcenter)+1];
    xe_interp(1,:,:,:) = importdata([pathname_Data1,'xe_',num2str(zint(1)),ID,'.mat']);
    xe_interp(2,:,:,:) = importdata([pathname_Data1,'xe_',num2str(zint(2)),ID,'.mat']);
    xe = min(max(xe_grid),max(min(xe_grid),squeeze(interp1(log(1+zint),xe_interp,log(1+zcenter)))));
    %%% END CHANGE
    
    for ii= 1:Nshells  %integral  
        R=Lpix*(Rmin(ii)+Rmax(ii))/2; % comoving radius of ring
        if R<70
            fcoll = LWgetGasShell2(fcoll_matrix,0,Rmax(ii))/(4*pi*Rmax(ii)^3/3);
            %%% MATAN CHANGE - 2015/10/06
            xeR0 = LWgetGasShell2(squeeze(fftn(xe)),0,Rmax(ii))./(4*pi*Rmax(ii)^3/3);        
            Maxfcoll = max(fcoll+xeR0/zeta,Maxfcoll);
            %Maxfcoll = max(fcoll,Maxfcoll); 
            %%% END CHANGE
        end
    end
    fcoll_matrix=[];

    PrevNeut=ones(128,128,128);
    if zcenter<40 && photoheatingVersion==2
        zs=[6:0.1:15 16:40];
        q=find(round(zs*10)/10==round(zcenter*10)/10);
        PrevNeut=logical(importdata([pathname_Data2,'xHI_',num2str(zs(q+1)),ID,'.mat'])>0);
    end
    
    %%% MATAN CHANGE - 2015/10/06
    Maxfcoll = max(fgas+xe/zeta,Maxfcoll);% adding the pixel
    %Maxfcoll = max(fgas,Maxfcoll);% adding the pixel
    Neut = (Maxfcoll<Threshold);% 1 if neutral, 0 if fully ionized
    xHI = max(0,(1-zeta*fgas-xe).*Neut.*PrevNeut);% ionized fraction of each pixel
    %xHI = max(0,(1-zeta*fgas).*Neut.*PrevNeut);% ionized fraction of each pixel
    %%% END CHANGE

    xHImean=mean(mean(mean(xHI)));
    
    save([pathname_Data2,'xHI_',num2str(zcenter),ID,'.mat'],'xHI');
    save([pathname_Data2,'Neut_',num2str(zcenter),ID,'.mat'],'Neut');
end