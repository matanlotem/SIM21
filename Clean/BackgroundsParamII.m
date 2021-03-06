%---calculates LyA, LW, Xrays radiative backgrounds, No LyA from Xrays at the moment
%---Tgas, xe, T-21

% Reion is the value of zeta

function [JLW21] = BackgroundsParamII(zcenter,ncube,fstar,flag,flagM,XeffTerm,Ispec,Reion,zeta,feedback,p,pop,FSfunc,photoheatingOn,photoheatingVersion)
%Reion = 0.075
 load Planck_parameters
%count = 0
%Lpix=3;
zMAX2 = 60;% start to calculate all radiative backgrounds and 21-cm
zMAX1 = 65;% start to evolve Tgas and xe
fstarM = fstar;% 0.1;
fstarA = fstar;%0.1;
global delta_cube
global pathname_Data1
global pathname_Data2
N=length(delta_cube);
JLW21 = zeros(N,N,N);
Lion = zeros(N,N,N);
eps = zeros(N,N,N);
Jalpha = zeros(N,N,N);


%save('test1z.mat','zMAX2');

if(zcenter>zMAX2)
   JLW21=1e-10*ones(N,N,N);
     save(strcat(pathname_Data1,'JLW_',num2str(zcenter),'_',num2str(ncube),'_',num2str(fstar),'_',num2str(flag),...
                    '_',num2str(flagM),'_',num2str(XeffTerm),'_',num2str(Ispec),'_',num2str(Reion),'_',num2str(feedback),...
                    '_',num2str(p),'_',num2str(pop),'_',num2str(FSfunc),'_',num2str(photoheatingVersion),'.mat'),'JLW21');
  
else       

%-----reionization parameters----------%
Threshold = 1/zeta;

%------------------------------------------%
   %save('test2z.mat','zMAX2');

    

    xHI_grid =1-[0,10^(-4),10^(-3.3), 10^(-3),10^(-2.6),10^(-2.3),10^(-2),10^(-1.6),10^(-1.3), 10^(-1),0.5,1];%2
    xe_grid =[10^(-4),10^(-3.3), 10^(-3),10^(-2.6),10^(-2.3),10^(-2)];% mean of the central pixel
    %zc_grid = 5:100;
    %Rset_grid = [linspace(1e-10,500/3,110) logspace(log10(550/3),log10(15000/3),20)];
    %R_grid = 3*(Rset_grid(1:end-1)+Rset_grid(2:end))/2;
             
IspecNew=Ispec;
    IsXRB=1;    %Q
    IsMQ=0;    %Q
    if (IspecNew>2) && (IspecNew<6)    %Q
        IsMQ=1;    %Q
        IspecNew=IspecNew-3;    %Q
    end    %Q
    if IspecNew==6    %Q
        IsXRB=0;    %Q
        IsMQ=1;     %Q 
        IspecNew=2;    %Q %Need to change the code - in this case we don't need to calculate the XRB's at all
    end    %Q
    
    OmegaMz=(1+(OLambda/Om)*(1+zcenter)^(-3))^(-1);
    d=OmegaMz-1;
    Dc=18*pi^2+82*d-39*d^2;
    MQzeta=Om/OmegaMz*Dc/(18*pi^2);
    
    C21 = 27*sqrt(0.15/Om/h^2)*(Ob*h^2/0.023);
    if(Lpix==3) 
       
        load('XCoefMatNEW_xHI0212.mat');%Q
        load('LionCoefMatNEW_xHI0212.mat');%Q
        XCoefMatMQ=XCoefMatNew/5*(MQzeta)^(5/6);%Q
        LionCoefMatMQ=LionCoefMatNew/5*(MQzeta)^(5/6);%Q
        XCoefMatNew=[];%Q
        LionCoefMatNew=[];%Q
        
       if(pop==2)
           load('LWCoefMatII.mat')
           LWCoefMat = LWCoefMatII;
           LWCoefMatII = []; 
           load('LyACoefMatII.mat')
           LyACoefMat = LyACoefMatII;
           LyACoefMatII = []; 
       else
           load('LWCoefMatIII.mat')
           LWCoefMat = LWCoefMatIII;
           LWCoefMatIII = []; 
           load('LyACoefMatIII.mat')
            LyACoefMat = LyACoefMatIII;
           LyACoefMatIII = []; 
       end
       
        if(IspecNew == 2)
            %----------------old power law SED--------------------------%
            %load('LionCoefMat_xHI2.mat');
            %load('XCoefMat_xHI2.mat'); 
            load('LionCoefMat_xHI0212.mat');
            load('XCoefMat_xHI0212.mat');
        else
            %----------------mixed X-ray SED--------------------------%
            
            %load('LionCoefMatNEW_xHI2.mat');% LionCoefMatNew; 
            %load('LionCoefMatNoA_xHI.mat');% LionCoefMatNew;
            load('LionCoefMatNEW_xHI0212.mat');% LionCoefMatNew; 
            load('LionCoefMatNoA_xHI0212.mat');% LionCoefMatNew;
            LionCoefMat = IspecNew*LionCoefMatNew+(1-IspecNew)*LionCoefMatNoA;
            LionCoefMatNew = [];
            LionCoefMatNoA = [];
            %load('XCoefMatNEW_xHI2.mat'); 
            %load('XCoefMatNoA_xHI.mat');l
            load('XCoefMatNEW_xHI0212.mat'); 
            load('XCoefMatNoA_xHI0212.mat');
            XCoefMat = IspecNew*XCoefMatNew+(1-IspecNew)*XCoefMatNoA;
            XCoefMatNew = []; 
            XCoefMatNoA = []; 
            % not updated
            load('LyAXCoefMat5000.mat');
            load('LyAXCoefMatNEW5000.mat'); 
            LyAXCoefMat = IspecNew*LyAXCoefMatNew+(1-IspecNew)*LyAXCoefMat;
            LyAXCoefMatNew = [];
            LyAXCoefMat = zeros(96,129);
        end
    else
        
        load('LWCoefMat6.mat');
        load('LyACoefMat6.mat');
        %----------------mixed X-ray SED--------------------------%
        load('LionCoefMatNEW6_xHI.mat');% LionCoefMatNew; 
        load('LionCoefMatNoA6_xHI.mat');% LionCoefMatNew;
        load('LionCoefMat6_xHI.mat');
       % LionCoefMat =  Ispec*LionCoefMatNew+(1-Ispec)*LionCoefMat; 
        LionCoefMat = IspecNew*LionCoefMatNew+(1-IspecNew)*LionCoefMatNoA;
        LionCoefMatNew = [];
        LionCoefMatNoA = [];
        load('XCoefMatNEW6_xHI.mat'); 
        load('XCoefMatNoA6_xHI.mat');
        load('XCoefMat6_xHI.mat');
        %XCoefMat = Ispec*XCoefMatNew+(1-Ispec)*XCoefMat;
        XCoefMat = IspecNew*XCoefMatNew+(1-IspecNew)*XCoefMatNoA;
        XCoefMatNew = []; 
        XCoefMatNoA = []; 
        % not updated
        % load('LyAXCoefMat5000.mat');
        % load('LyAXCoefMatNEW5000.mat'); 
        LyAXCoefMat = zeros(96,129);

    end

 %save('test3z.mat','zMAX2');
    zc = 6:100;
    z=(6:1:75);
    dz=0.0001;
    zp = z+dz;

    Rset = [linspace(1e-10,500/3,110) logspace(log10(550/3),log10(15000/3),20)];
    Rmin = Rset(1:end-1);
    Rmax = Rset(2:end);
    Nshells = length(Rmin);
    if(length(Rmax) ~= length(Rmin))
        pause;
    end

%if (zcenter < zMAX2+1)
    zsM =max(max(getzmax(zcenter,2 ),getRtoz(140,zcenter)),65);
    Ind = find(z==zcenter): find(z==ceil(zsM));
    dfdt_matrix = zeros(length(Ind),N,N,N); 
    if Ispec>2
    dfdtMQ_matrix = zeros(length(Ind),N,N,N); 
    end
   % fcoll_matrix = zeros(length(Ind),N,N,N);  
    fcoll_matrix = zeros(N,N,N); 
    xe_matrix = zeros(N,N,N); 
    xHI_matrixCurrent=zeros(N,N,N);
    xHI_matrixTemp = zeros(3,N,N,N); 
   % Neut_matrix = zeros(length(Ind),N,N,N); 
    for ii= 1:length(Ind)
        zii = z(Ind(ii));
        D = LWgetDz(zii)/D40;
        Dp = LWgetDz(zii+dz)/D40;    
        if(feedback)
            z0 = (1+zii)/(min(max(p,0.05),1))^(2/3)-1;
            J_interp = zeros(2,N,N,N);
            zint = (ceil(z0)-1:ceil(z0));

            % load the LW to calculate the feedback
            for jj=1:2    
                
                load(strcat(pathname_Data1,'JLW_',num2str(zint(jj)),'_',num2str(ncube),'_',num2str(fstar),'_',num2str(flag),...
                    '_',num2str(flagM),'_',num2str(XeffTerm),'_',num2str(Ispec),'_',num2str(Reion),'_',num2str(feedback),...
                    '_',num2str(p),'_',num2str(pop),'_',num2str(FSfunc),'_',num2str(photoheatingVersion),'.mat'));% total LyA flux
                J_interp(jj,:,:,:) = JLW21; 
                JLW21 = [];
            end      
            JLW21 = squeeze(exp(interp1(zint,log(abs(J_interp)),z0))); 
           
        end
        [fgas_z,fgas_zp,fgasMQ_z,fgasMQ_zp] = grid_interpSF2(flag,flagM,feedback*JLW21,zii,Ispec,fstarM,fstarA,FSfunc,photoheatingOn,photoheatingVersion,1);% fstar inside %Q
                
        dt = getztot(zii+dz)-getztot(zii);% [years]    
        dfdt_matrix(ii,:,:,:) = fftn((fgas_zp.*(1+delta_cube*Dp)-fgas_z.*(1+delta_cube*D)).*Lpix^3/dt);% 
        if Ispec>2
        dfdtMQ_matrix(ii,:,:,:) = fftn((fgasMQ_zp.*(1+delta_cube*Dp)-fgasMQ_z.*(1+delta_cube*D)).*Lpix^3/dt);%  %Q
        end
        
        %fcoll_matrix(ii,:,:,:) =fftn(fgas_z/fstar);
        if(ii==1)
             fcoll_matrix =fftn(fgas_z/fstar);
        end
        fgas_z = [];
        fgas_zp = [];
%         if(zii>zMAX2)
%             xe =  (1.716e-11*zii.^3-3.823e-9*zii.^2+4.941e-7*zii+0.0001056).*ones(N,N,N);
%             xe_matrix(ii,:,:,:) =fftn(xe);
%         else
%             load(strcat(pathname_Data1,'xe_',num2str(zii),'_',num2str(ncube),'_',num2str(fstar),'_',num2str(flag),...
%                     '_',num2str(flagM),'_',num2str(XeffTerm),'_',num2str(Ispec),'_',num2str(Reion),'_',num2str(feedback),...
%                     '_',num2str(p),'_',num2str(pop),'.mat'));% total LyA flux
%             xe_matrix(ii,:,:,:) = fftn(xe);%xe at each pixel at zii
%             xe = [];
%         end
    end

     for ii= 1:3
        zii = z(Ind(ii));
        if(zii>zcenter)
            if(zii>zMAX2)
                xHI_matrixTemp(ii,:,:,:) = fftn(ones(N,N,N));
                %Neut_matrix(ii,:,:,:) =fftn( ones(N,N,N));
            else
                load(strcat(pathname_Data2,'xHI_',num2str(zii),'_',num2str(ncube),'_',num2str(fstar),'_',num2str(flag),...
                    '_',num2str(flagM),'_',num2str(XeffTerm),'_',num2str(Ispec),'_',num2str(Reion),'_',num2str(feedback),'_',...
                    num2str(p),'_',num2str(pop),'_',num2str(FSfunc),'_',num2str(photoheatingVersion),'.mat'));% total LyA flux
          
                xHI_matrixTemp(ii,:,:,:) =fftn(max(0,min(1,xHI)));%xe at each pixel at zii
                xHI = [];
            end
        end
     end
    
         
         
    zextrap = z(Ind);
    xHI_matrixCurrent(:,:,:) =squeeze(interp1(log(1+zextrap(2:3)),xHI_matrixTemp(2:3,:,:,:),log(1+zextrap(1)),'linear','extrap'));
    xHI_matrixTemp=[];
    
    %count = 1
    %save('test4z.mat','zMAX2');
    dfdt_interp = zeros(2,N,N,N);
    xe_interp = zeros(2,N,N,N);  
    xHI_interp = zeros(2,N,N,N);
    Maxfcoll = zeros(N,N,N);
    eps = zeros(N,N,N);
    Lion = zeros(N,N,N);
    Jalpha=zeros(N,N,N);
    %JAX=zeros(N,N,N);
    JA=zeros(N,N,N);
    JLW21=zeros(N,N,N);
    for ii= 1:Nshells  %integral
   
        R=Lpix*(Rmin(ii)+Rmax(ii))/2; % comoving radius of ring
        zshell = getRtoz(R,zcenter);
        if(zshell<ceil(zsM)&(zshell>zcenter))
            izs = find(z>zshell,1,'first');
            i_interp=0;
            for iz=izs-1:izs
                i_interp=i_interp+1;
                indZ=find(Ind==find(z==z(iz)));
                dfdt_interp(i_interp,:,:,:) = LWgetGasShell2(squeeze(dfdt_matrix(indZ,:,:,:)),Rmin(ii),Rmax(ii));% df/dt in [sec^-1] units
                if Ispec>2
                dfdtMQ_interp(i_interp,:,:,:) = LWgetGasShell2(squeeze(dfdtMQ_matrix(indZ,:,:,:)),Rmin(ii),Rmax(ii));% df/dt in [sec^-1] units %Q
                end
               
                    if(z(iz)>zMAX2)
                                 xe =  (1.716e-11*z(iz).^3-3.823e-9*z(iz).^2+4.941e-7*z(iz)+0.0001056).*ones(N,N,N);
                                 xe_matrix(:,:,:) =fftn(xe);
                          else
                                 load(strcat(pathname_Data1,'xe_',num2str(z(iz)),'_',num2str(ncube),'_',num2str(fstar),'_',num2str(flag),...
                                '_',num2str(flagM),'_',num2str(XeffTerm),'_',num2str(Ispec),'_',num2str(Reion),'_',num2str(feedback),...
                                '_',num2str(p),'_',num2str(pop),'_',num2str(FSfunc),'_',num2str(photoheatingVersion),'.mat'));% total LyA flux
                                xe_matrix(:,:,:) = fftn(xe);%xe at each pixel at zii
                                xe = [];
                          end
                
                 xe_interp(i_interp,:,:,:) = LWgetGasShell2(squeeze(xe_matrix(:,:,:)),0,Rmax(ii))/(4*pi*Rmax(ii)^3/3); 
                 
                 
               if(z(iz)>zcenter)
                            if(z(iz)>zMAX2)
                                xHI_matrix(:,:,:) = fftn(ones(N,N,N));
                                %Neut_matrix(ii,:,:,:) =fftn( ones(N,N,N));
                            else
                                load(strcat(pathname_Data2,'xHI_',num2str(z(iz)),'_',num2str(ncube),'_',num2str(fstar),'_',num2str(flag),...
                                    '_',num2str(flagM),'_',num2str(XeffTerm),'_',num2str(Ispec),'_',num2str(Reion),'_',num2str(feedback),'_',...
                                    num2str(p),'_',num2str(pop),'_',num2str(FSfunc),'_',num2str(photoheatingVersion),'.mat'));% total LyA flux

                                xHI_matrix(:,:,:) =fftn(max(0,min(1,xHI)));%xe at each pixel at zii
                                xHI = [];
                            end
                        end

                        if z(iz)==zcenter
                            xHI_matrix=xHI_matrixCurrent;
                        end

                
                xHI_interp(i_interp,:,:,:) = LWgetGasShell2(squeeze(xHI_matrix(:,:,:)),0,Rmax(ii))./(4*pi*Rmax(ii)^3/3); 
           
            end
            
            aa=1;
            xeshell =min(max(xe_grid),max(min(xe_grid),squeeze(interp1(log(1+z(izs-1:izs)),xe_interp,log(1+zshell)))));
            %mean(mean(mean(xeshell)))
            xHIshell = min(max(xHI_grid),max(min(xHI_grid),squeeze(interp1(log(1+z(izs-1:izs)),xHI_interp,log(1+zshell)))));
            %mean(mean(mean(xHIshell)))
            dfdt =squeeze(exp(interp1(log(1+z(izs-1:izs)),log(abs(dfdt_interp)+1e-26),log(1+zshell))));   %dfdt per shell (including the shell volume)
            dfdtMQ=zeros(size(dfdt));
            if Ispec>2
            dfdtMQ =squeeze(exp(interp1(log(1+z(izs-1:izs)),log(abs(dfdtMQ_interp)),log(1+zshell))));   %Q  
            end
            
           
            
            JA =  JA+dfdt*LyACoefMat(find(zc==zcenter),ii);% 
            eps =  eps+IsXRB.*XeffTerm*dfdt.*10.^(interp2(log10(xHI_grid+1e-16),log10(xe_grid+1e-16),...
               log10(squeeze(XCoefMat(find(zc==zcenter),ii,:,:))),log10(squeeze(xHIshell)+1e-16),...
               log10(squeeze(xeshell)+1e-16))) + IsMQ.*dfdtMQ.*XeffTerm.*(0.05./fstar).*((1+zshell)/10).*10.^(interp2(log10(xHI_grid+1e-16),log10(xe_grid+1e-16),...
                log10(squeeze(XCoefMatMQ(find(zc==zcenter),ii,:,:))),log10(squeeze(xHIshell)+1e-16),...
                log10(squeeze(xeshell)+1e-16)));% eV/sec/baryon  %Q
            
          
            Lion =  Lion+IsXRB.*XeffTerm*dfdt.*10.^(interp2(log10(xHI_grid+1e-16),log10(xe_grid+1e-16),...
                log10(squeeze(LionCoefMat(find(zc==zcenter),ii,:,:))),log10(squeeze(xHIshell)+1e-16),...
                log10(squeeze(xeshell)+1e-16))) + IsMQ.*dfdtMQ.*XeffTerm.*(0.05./fstar).*((1+zshell)/10).*10.^(interp2(log10(xHI_grid+1e-16),log10(xe_grid+1e-16),...
                log10(squeeze(LionCoefMatMQ(find(zc==zcenter),ii,:,:))),log10(squeeze(xHIshell)+1e-16),...
                log10(squeeze(xeshell)+1e-16)));% eV/sec/baryon  %Q

            
%           JAX = JAX+ XeffTerm*dfdt.*10.^(interp2(log10(xHI_grid+1e-16),log10(xe_grid+1e-16),...
%                 log10(squeeze(LyAXCoefMat(find(zc==zcenter),ii,:,:))),log10(squeeze(xHIshell)+1e-16),...
%                 log10(squeeze(xeshell)+1e-16)));
            
          JLW21 =  JLW21+dfdt*LWCoefMat(find(zc==zcenter),ii);% 
     
            
            %-------------reionization------------------% 
           if(R<70) 
                fcoll = LWgetGasShell2(fcoll_matrix,0,Rmax(ii))/(4*pi*Rmax(ii)^3/3);
                Maxfcoll = max(fcoll,Maxfcoll); 
%             mean(mean(mean(Maxfcoll)))
%             R
           end
           dfdt = [];  
           fcoll = []; 
        end
    end
   %save('test5z.mat','zMAX2');
   dfdt_matrix=[];
   dfdt_interp = [];
   fcoll_matrix=[];
   fcoll_interp = [];
   xe_interp=[];
   xe_matrix=[];
   xHI_matrix=[];
   xHI_interp=[];
   idfdtshell=[];
   
   %JAX=  zeros(N,N,N);%JAX.*(1+delta_cube.*LWgetDz(zcenter)/LWgetDz(40)).*nb*(1+zcenter)^3*c/(getHz(zcenter)*4*pi*nu_a^2*hpl)/(3*10^24)^2;
   Jalpha = JA;%+JAX;
   JLW21=JLW21./10^(-21);
   
    save(strcat(pathname_Data1,'Lion_',num2str(zcenter),'_',num2str(ncube),'_',num2str(fstar),'_',num2str(flag),...
                        '_',num2str(flagM),'_',num2str(XeffTerm),'_',num2str(Ispec),'_',num2str(Reion),'_',num2str(feedback),'_',...
                       num2str(p),'_',num2str(pop),'_',num2str(FSfunc),'_',num2str(photoheatingVersion),'.mat'),'Lion');
    save(strcat(pathname_Data1,'eps_',num2str(zcenter),'_',num2str(ncube),'_',num2str(fstar),'_',num2str(flag),...
                        '_',num2str(flagM),'_',num2str(XeffTerm),'_',num2str(Ispec),'_',num2str(Reion),'_',num2str(feedback),...
                        '_',num2str(p),'_',num2str(pop),'_',num2str(FSfunc),'_',num2str(photoheatingVersion),'.mat'),'eps');
    save(strcat(pathname_Data1,'Jalpha_',num2str(zcenter),'_',num2str(ncube),'_',num2str(fstar),'_',num2str(flag),...
                        '_',num2str(flagM),'_',num2str(XeffTerm),'_',num2str(Ispec),'_',num2str(Reion),'_',num2str(feedback),...
                        '_',num2str(p),'_',num2str(pop),'_',num2str(FSfunc),'_',num2str(photoheatingVersion),'.mat'),'Jalpha');
    save(strcat(pathname_Data1,'JLW_',num2str(zcenter),'_',num2str(ncube),'_',num2str(fstar),'_',num2str(flag),...
                         '_',num2str(flagM),'_',num2str(XeffTerm),'_',num2str(Ispec),'_',num2str(Reion),'_',num2str(feedback),'_',...
                         num2str(p),'_',num2str(pop),'_',num2str(FSfunc),'_',num2str(photoheatingVersion),'.mat'),'JLW21');

                
   Jalpha = [];
   JA = [];
   %JAX = [];
  JLW21 = [];
   
   %%% JLW to use as input for ionization fraction calculation ------
   
    if (feedback==1)
        zii =  zcenter;
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
   %%%-----------------------------
    
   PrevNeut=ones(N,N,N);
   if zcenter<40 && photoheatingVersion==2
   zs=[6:0.1:15 16:40];
   q=find(round(zs*10)/10==round(zcenter*10)/10);
   load(strcat(pathname_Data2,'xHI_',num2str(zs(q+1)),'_',num2str(ncube),'_',num2str(fstar),'_',num2str(flag),...
                        '_',num2str(flagM),'_',num2str(XeffTerm),'_',num2str(Ispec),'_',num2str(Reion),'_',num2str(feedback),'_',...
                        num2str(p),'_',num2str(pop),'_',num2str(FSfunc),'_',num2str(photoheatingVersion),'.mat'),'xHI');
    PrevNeut=logical(xHI>0);
    end 
   
   
    fgas = grid_interpSF2(flag,flagM,feedback*JLW21,zcenter,Ispec,fstar,fstar,FSfunc,photoheatingOn,photoheatingVersion,0)/fstar;
    Maxfcoll = max(fgas,Maxfcoll);% adding the pixel
    Neut = (Maxfcoll<Threshold);% 1 if neutral, 0 if fully ionized
    xHI = max(0, (1-zeta*fgas).*Neut.*PrevNeut);% neutral fraction of each pixel (as in 21cmFAST)

     save(strcat(pathname_Data2,'xHI_',num2str(zcenter),'_',num2str(ncube),'_',num2str(fstar),'_',num2str(flag),...
                        '_',num2str(flagM),'_',num2str(XeffTerm),'_',num2str(Ispec),'_',num2str(Reion),'_',num2str(feedback),...
                        '_',num2str(p),'_',num2str(pop),'_',num2str(FSfunc),'_',num2str(photoheatingVersion),'.mat'),'xHI');

     save(strcat(pathname_Data2,'Neut_',num2str(zcenter),'_',num2str(ncube),'_',num2str(fstar),'_',num2str(flag),...
                    '_',num2str(flagM),'_',num2str(XeffTerm),'_',num2str(Ispec),'_',num2str(Reion),'_',num2str(feedback),...
                    '_',num2str(p),'_',num2str(pop),'_',num2str(FSfunc),'_',num2str(photoheatingVersion),'.mat'),'Neut');
          
% ------------- find Tspin and 21-cm ------------------%    
%------------------------------
    load(strcat(pathname_Data1,'xe_',num2str(zcenter),'_',num2str(ncube),'_',num2str(fstar),'_',num2str(flag),...
                        '_',num2str(flagM),'_',num2str(XeffTerm),'_',num2str(Ispec),'_',num2str(Reion),'_',num2str(feedback),...
                        '_',num2str(p),'_',num2str(pop),'_',num2str(FSfunc),'_',num2str(photoheatingVersion),'.mat'));

    load(strcat(pathname_Data2,'TK_',num2str(zcenter),'_',num2str(ncube),'_',num2str(fstar),'_',num2str(flag),...
                        '_',num2str(flagM),'_',num2str(XeffTerm),'_',num2str(Ispec),'_',num2str(Reion),'_',num2str(feedback),...
                        '_',num2str(p),'_',num2str(pop),'_',num2str(FSfunc),'_',num2str(photoheatingVersion),'.mat'));
% count = 2  

            %save('test6z.mat','zMAX2');
    Ts = getTs(TK, zcenter,ncube,fstar,flag,flagM,XeffTerm,Ispec,Reion,feedback,p,pop,FSfunc,photoheatingVersion);
    Ts(Ts==0)=1e-20;   
    deltaTerm = (ones(size(delta_cube))+max(-0.9,min(1,delta_cube*LWgetDz(zcenter)/LWgetDz(40))));
    T21cm = C21*sqrt((1+zcenter)/10).*deltaTerm.*(1-2.725*(1+zcenter)./Ts).*xHI;
     save(strcat(pathname_Data2,'T21cm_',num2str(zcenter),'_',num2str(ncube),'_',num2str(fstar),'_',num2str(flag),...
                        '_',num2str(flagM),'_',num2str(XeffTerm),'_',num2str(Ispec),'_',num2str(Reion),'_',num2str(feedback),...
                        '_',num2str(p),'_',num2str(pop),'_',num2str(FSfunc),'_',num2str(photoheatingVersion),'.mat'),'T21cm');

    T21cm=[];
    Ts=[];
    deltaTerm=[];
end
 %count = 3

%-----------------------------------------------------------------------%
%------------ gas temperature and xe -----------------------------------%
if (zcenter < zMAX1+1)
    if (zcenter > zMAX2)
         TK = (1+zcenter).^2*41.6./(1+45).^2.*ones(N,N,N); 
         xe = (1.716e-11*zcenter.^3-3.823e-9*zcenter.^2+4.941e-7*zcenter+0.0001056).*ones(N,N,N);
    end

%save('test7z.mat','zMAX2');
%--------------Runge-Kutta-----------------------------------------%

%---- step1
    x2 = log(1+zcenter-1);
    x0 = log(1+zcenter);
    A0 = log(1+xe);% at zcenter
    B0 = log(TK);
    xe=[];
    TK=[];
    Lion0=Lion;
    eps0=eps;

%--- step2
    x1 = (x0+x2)/2;
    z1 = exp(x1)-1;
    A1 = A0+getf(zcenter,A0,B0,Lion0).*(x1-x0);
    B1 = B0+getg(zcenter,A0,B0,Lion0,eps0).*(x1-x0);

    zi = [zcenter,zcenter+1];
    Lion_extrap=zeros(2,N,N,N);
    eps_extrap=zeros(2,N,N,N);
    for indz=1:length(zi)
        if (zi(indz) == zcenter)
            Lion_extrap(indz,:,:,:) = Lion;  
            eps_extrap(indz,:,:,:) = eps;  
        elseif (zi(indz)>zMAX2) 
            Lion_extrap(indz,:,:,:)=zeros(N,N,N);
            eps_extrap(indz,:,:,:)=zeros(N,N,N);
        else
            load(strcat(pathname_Data1,'Lion_',num2str(zi(indz)),'_',num2str(ncube),'_',num2str(fstar),'_',num2str(flag),...
                        '_',num2str(flagM),'_',num2str(XeffTerm),'_',num2str(Ispec),'_',num2str(Reion),'_',num2str(feedback),...
                        '_',num2str(p),'_',num2str(pop),'_',num2str(FSfunc),'_',num2str(photoheatingVersion),'.mat'),'Lion');
            load(strcat(pathname_Data1,'eps_',num2str(zi(indz)),'_',num2str(ncube),'_',num2str(fstar),'_',num2str(flag),...
                        '_',num2str(flagM),'_',num2str(XeffTerm),'_',num2str(Ispec),'_',num2str(Reion),'_',num2str(feedback),...
                        '_',num2str(p),'_',num2str(pop),'_',num2str(FSfunc),'_',num2str(photoheatingVersion),'.mat'),'eps');

            Lion_extrap(indz,:,:,:) = Lion;% df/dt in [sec^-1] units 
           eps_extrap(indz,:,:,:) = eps;% df/dt in [sec^-1] units    
        end
    end
    Lion1 = squeeze(interp1(log10(1+zi), Lion_extrap,log10(1+z1),'linear','extrap'));% extrap Lion at the intermediate step
    eps1 =  squeeze(interp1(log10(1+zi), eps_extrap,log10(1+z1),'linear','extrap'));% extrap Lion at the intermediate step
    Lion_extrap = [];
    eps_extrap=[];
%-----step 3

    A2 = A0+ getf(z1,A1,B1,Lion1).*(x2-x0);
    B2 = B0+ getg(z1,A1,B1,Lion1,eps1).*(x2-x0);

    xe = exp(A2)-1;
    TK = exp(B2);
     save(strcat(pathname_Data1,'xe_',num2str(zcenter-1),'_',num2str(ncube),'_',num2str(fstar),'_',num2str(flag),...
                        '_',num2str(flagM),'_',num2str(XeffTerm),'_',num2str(Ispec),'_',num2str(Reion),'_',num2str(feedback),...
                        '_',num2str(p),'_',num2str(pop),'_',num2str(FSfunc),'_',num2str(photoheatingVersion),'.mat'),'xe');
     save(strcat(pathname_Data2,'TK_',num2str(zcenter-1),'_',num2str(ncube),'_',num2str(fstar),'_',num2str(flag),...
                        '_',num2str(flagM),'_',num2str(XeffTerm),'_',num2str(Ispec),'_',num2str(Reion),'_',num2str(feedback),...
                        '_',num2str(p),'_',num2str(pop),'_',num2str(FSfunc),'_',num2str(photoheatingVersion),'.mat'),'TK');

%


  %  count = 4                  
TK = [];
xe = [];
Lion1=[];
eps1=[];
Lion = [];
eps = [];
end
end