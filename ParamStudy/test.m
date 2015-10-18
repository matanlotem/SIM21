addpath('/a/home/cc/tree/taucc/students/physics/matanlotem/Work/SIM21/lib/');
%p = paramStudy();
%pt = paramStudy();
%for caseNum = [1:5]
%	c = pt.paramCases(caseNum).c;
%	runName = ['PS_New_',num2str(caseNum)];
%	c.setPathExt([runName,'_New']);
%	c.setOutputPath([c.outputPath(1:end-1),'_New/'],1)
%end

%for caseNum = [1:5]
%	c = p.paramCases(caseNum).c;
%	ct = pt.paramCases(caseNum).c;
%	SIM21Analysis.plotAllGraphsByZ([c,ct]);
%end
%
%pt.plot4CasesZGraphs();



tc1 = SIM21Case('MLrun1','AFtest1','/scratch300/matanlotem/Temp/Output/',21,0.05,1,35.5,588.2,2,42.5,0,0,2,1,2,42.5);
tc2 = SIM21Case('MLrun2','AFtest2','/scratch300/matanlotem/Temp/Output/',21,0.05,1,16.5,1258.3,2,72,0,0,2,1,2,72);
tc3 = SIM21Case('MLrun3','AFtest3','/scratch300/matanlotem/Temp/Output/',21,0.05,1,35.5,682.8,6,180,0,0,2,1,2,180);

tc1b = SIM21Case('MLrunZfix1','AFtest1b','/scratch300/matanlotem/Temp/Output_b/',21,0.05,1,35.5,588.2,2,42.5,0,0,2,1,2,42.5);
tc2b = SIM21Case('MLrunZfix2','AFtest2b','/scratch300/matanlotem/Temp/Output_b/',21,0.05,1,16.5,1258.3,2,72,0,0,2,1,2,72);
tc3b = SIM21Case('MLrunZfix3','AFtest3b','/scratch300/matanlotem/Temp/Output_b/',21,0.05,1,35.5,682.8,6,180,0,0,2,1,2,180);

tc3c = SIM21Case('MLrunMQfix3','AFtest3c','/scratch300/matanlotem/Temp/Output_c/',21,0.05,1,35.5,682.8,6,180,0,0,2,1,2,180);


ac1 = copy(tc1);
ac1.name = 'AFrun1';
ac1.dataPath = '/scratch300/Anastasia_F/Run_DATA/';
ac1.tmpDataPath = '/scratch/prun/P_DATA/'; % compute-0-65
ac1.setOutputPath('/scratch300/matanlotem/Temp/Output_a/')
ac2 = copy(tc2);
ac2.name = 'AFrun2';
ac2.dataPath = '/scratch300/Anastasia_F/Run_DATA/';
ac2.setOutputPath('/scratch300/matanlotem/Temp/Output_a/')
ac3 = copy(tc3);
ac3.name = 'AFrun3';
ac3.dataPath = '/scratch300/Anastasia_F/Run_DATA/';
ac3.tmpDataPath = '/scratch/prun/P_DATA/'; % compute-0-65
ac3.setOutputPath('/scratch300/matanlotem/Temp/Output_a/')

SIM21Analysis.plotAllGraphsByZ([ac3,tc3,tc3b,tc3c]);


%%tc1 = SIM21Case('AFtest1b','AFtest1b','/scratch300/matanlotem/Temp/Output/',21,0.05,1,35.5,588.2,2,42.5,0,0,2,1,2,42.5);
%%tc2 = SIM21Case('AFtest2b','AFtest2b','/scratch300/matanlotem/Temp/Output/',21,0.05,1,16.5,1258.3,2,72,0,0,2,1,2,72);
%%tc3 = SIM21Case('AFtest3b','AFtest3b','/scratch300/matanlotem/Temp/Output/',21,0.05,1,35.5,682.8,6,180,0,0,2,1,2,180);
%%
%%
%%ac1 = copy(tc1);
%%ac1.dataPath = '/scratch300/Anastasia_F/Run_DATA/';
%%ac1.tmpDataPath = '/scratch/prun/P_DATA/'; % compute-0-65
%%ac1.setOutputPath('/scratch300/matanlotem/Temp/Output_a/')
%%ac2 = copy(tc2);
%%ac2.dataPath = '/scratch300/Anastasia_F/Run_DATA/';
%%ac2.setOutputPath('/scratch300/matanlotem/Temp/Output_a/')
%%ac3 = copy(tc3);
%%ac3.dataPath = '/scratch300/Anastasia_F/Run_DATA/';
%%ac3.tmpDataPath = '/scratch/prun/P_DATA/'; % compute-0-65
%%ac3.setOutputPath('/scratch300/matanlotem/Temp/Output_a/')
%%
%%
%%cd('~/../anstasi3/work/BackgroundsCODE_BigBox_2/');
%%global pathname;
%%global pathname_Data1;
%%global pathname_Data2;
%%global pathname_Data3;
%%global pathname_Code;
%%
%%pathname = '/home/fialkov/work/';
%%pathname_Data1 = '/scratch/prun/P_DATA/';
%%pathname_Data2 = '/scratch300/Anastasia_F/Run_DATA/';
%%pathname_Data3 = '/scratch300/Anastasia_F/';
%%pathname_Code = '/home/fialkov/work/21cmCODE/';
%%
%%ncube = ac1.ncube;
%%fstar = ac1.fstar;
%%flag = ac1.vbc;
%%flagM = ac1.vc;
%%XeffTerm = ac1.fx;
%%Ispec = ac1.sed;
%%Reion = ac1.tau;
%%zeta = ac1.zeta;
%%feedback = ac1.feedback;
%%p = ac1.delayParam;
%%pop = ac1.pop;
%%FSfunc = ac1.fsfunc;
%%photoheatingVersion = ac1.phVersion;
%%
%%
%%global delta_cube;
%%global vbc_cube;
%%global gridInterpKey;
%%delta_cube = importdata(strcat(pathname_Data3,'DataBackgrounds_withPlanck/my',num2str(ncube),'_d.dat'));
%%vbc_cube = importdata(strcat(pathname_Data3,'DataBackgrounds_withPlanck/my',num2str(ncube),'_v.dat'));
%%gridInterpKey = ac1.ID;
%%
%%%zcenter = 60;
%%zcenter = 14.9;
%%photoheatingOn = photoheatingVersion ~= 0;
%%load Planck_parameters;
%%
%%xe_grid =[10^(-4),10^(-3.3), 10^(-3),10^(-2.6),10^(-2.3),10^(-2),10^(-1.6),10^(-1.3), 10^(-1),0.5,0.9,0.99,1];% mean of the central pixel
%%Rset = [linspace(1e-10,500/3,110) logspace(log10(550/3),log10(15000/3),20)];
%%Rmin = Rset(1:end-1);
%%Rmax = Rset(2:end);
%%Nshells = length(Rmin);
%%if(length(Rmax) ~= length(Rmin))
%%    pause;
%%end
%%
%%%zfgas = 5:100;
%%zgas = (6:1:16);
%%JLW21 = zeros(N,N,N);
%%
%%fgas = grid_interpSFAviad2(flag,flagM,feedback*JLW21,zcenter,fstar,fstar,zeta,FSfunc,photoheatingOn,photoheatingVersion);
%%fgas=fgas/fstar;
%%fcoll_matrix =fftn(fgas);
%%Maxfcoll = zeros(N,N,N);
%%
%%%___ 22.09.2015 AF added
%%% added to do x_e properly
%%zint = [floor(zcenter),floor(zcenter)+1];
%%xe_interp = zeros(2,N,N,N);
%%for kk = 1:2
%%                load(strcat(pathname_Data1,'xe_',num2str(zint(kk)),'_',num2str(ncube),'_',num2str(fstar),'_',num2str(flag),...
%%                                        '_',num2str(flagM),'_',num2str(XeffTerm),'_',num2str(Ispec),'_',num2str(Reion),'_',num2str(feedback),...
%%                                        '_',num2str(p),'_',num2str(pop),'_',num2str(FSfunc),'_',num2str(photoheatingVersion),'.mat'));
%%                xe_interp(kk,:,:,:) = xe; 
%%                xe =[]; 
%%end
%%xepix0  =min(max(xe_grid),max(min(xe_grid),squeeze(interp1(log(1+zint),xe_interp,log(1+zcenter)))));
%%xepix = fftn(xepix0);
%%               
%xeR = zeros(N,N,N);
%%for ii= 1:Nshells  %integral
%%        R=Lpix*(Rmin(ii)+Rmax(ii))/2; 
%%        
%%         if(R<70) 
%%                    fcoll = LWgetGasShell2(fcoll_matrix,0,Rmax(ii))/(4*pi*Rmax(ii)^3/3);
%%                    xeR0 = LWgetGasShell2(squeeze(xepix),0,Rmax(ii))./(4*pi*Rmax(ii)^3/3);        
%%                    Maxfcoll = max(fcoll+xeR0/zeta,Maxfcoll);
%%                 %   IndR = find(Maxfcoll==fcoll);
%%               %     xeR(IndR) = xepix(IndR);% mean fraction of ionized within radius R
%%         end          
%%end
%%%%%%
%%%___
%%
%%% load the result from zcenter+0.1
%%PrevNeut=ones(N,N,N);
%%if zcenter<40 && photoheatingVersion==2
%%   zs=[6:0.1:15 16:40];
%%   q=find(round(zs*10)/10==round(zcenter*10)/10);
%%   load(strcat(pathname_Data2,'xHI_',num2str(zs(q+1)),'_',num2str(ncube),'_',num2str(fstar),'_',num2str(flag),...
%%                        '_',num2str(flagM),'_',num2str(XeffTerm),'_',num2str(Ispec),'_',num2str(Reion),'_',num2str(feedback),'_',...
%%                        num2str(p),'_',num2str(pop),'_',num2str(FSfunc),'_',num2str(photoheatingVersion),'.mat'),'xHI');
%%    PrevNeut=logical(xHI>0);
%%end
%%
%%%___ 22.09.2015 AF added
%%
%%    
%%Maxfcoll = max(fgas+xepix0/zeta,Maxfcoll);% adding the pixel
%%Neut = (Maxfcoll<Threshold);% 1 if neutral, 0 if fully ionized
%%
%%%IndR = find(Maxfcoll==fgas);
%%%xeR(IndR) = xepix(IndR);% mean fraction of ionized within radius R
%%xHI = max(0, (1-zeta*fgas-xepix0).*Neut.*PrevNeut);% ionized fraction of each pixel
%%%xHI = max(0, (1-zeta*fgas).*Neut.*PrevNeut);% ionized fraction of each pixel
%%%___
%%
%%
%%xHImean=mean(mean(mean(xHI)));




%load('LyACoefMatII.mat')
%LyACoefMat = LyACoefMatII;
%LyACoefMatII = []; 
%load(strcat('LyAXCoefMat_xHI_160915.mat'),'LyAXCoefMat');
%load(strcat('XCoefMat_xHI_160915.mat'),'XCoefMat');
%load(strcat('LionCoefMat_xHI_160915.mat'),'LionCoefMat');       
%
%zMAX2 = 60;% start to calculate all radiative backgrounds and 21-cm
%zMAX1 = 65;% start to evolve Tgas and xe
%fstarM = fstar;% 0.1;
%fstarA = fstar;%0.1;
%N=length(delta_cube);
%JLW21 = zeros(N,N,N);
%Lion = zeros(N,N,N);
%eps = zeros(N,N,N);
%JAX = zeros(N,N,N);
%Jalpha = zeros(N,N,N);
%LionMQ = zeros(N,N,N);
%epsMQ = zeros(N,N,N);
%JAXMQ = zeros(N,N,N);
%
%MQ = Ispec > 3;
%Threshold = 1/zeta;
%xHI_grid =1-[0,10^(-4),10^(-3.3), 10^(-3),10^(-2.6),10^(-2.3),10^(-2),10^(-1.6),10^(-1.3), 10^(-1),0.5,1];%2
%xe_grid =[10^(-4),10^(-3.3), 10^(-3),10^(-2.6),10^(-2.3),10^(-2),10^(-1.6),10^(-1.3), 10^(-1),0.5,0.9,0.99,1];% mean of the central pixel
%
%C21 = 27*sqrt(0.15/Om/h^2)*(Ob*h^2/0.023);
%
%zc = 5:100;
%z=(5:1:75);
%dz=0.0001;
%zp = z+dz;
%Rset = [linspace(1e-10,500/3,110) logspace(log10(550/3),log10(15000/3),20)];
%Rmin = Rset(1:end-1);
%Rmax = Rset(2:end);
%Nshells = length(Rmin);
%
%zsM =max(max(getzmax(zcenter,2),getRtoz(140,zcenter)),65);
%Ind = find(z==zcenter): find(z==ceil(zsM));
%dfdt_matrix = zeros(length(Ind),N,N,N);
%for ii= 1:length(Ind)      
%    zii = z(Ind(ii));
%    D = LWgetDz(zii)/D40;
%    Dp = LWgetDz(zii+dz)/D40;
%	[fgas_z,fgas_zp,fgasMQ_z,fgasMQ_zp] = grid_interpSFAviad(flag,flagM,feedback*JLW21,zii,Ispec,fstarM,fstarA,zeta,FSfunc,photoheatingOn,photoheatingVersion,1);% fstar inside %Q
%	dt = getztot(zii+dz)-getztot(zii);% [years]    
%    dfdt_matrix(ii,:,:,:) = fftn((fgas_zp.*(1+delta_cube*Dp)-fgas_z.*(1+delta_cube*D)).*Lpix^3/dt);% 
%
%    if(ii==1)
%             fcoll_matrix =fftn(fgas_z/fstar);
%	end
%end
%
%for ii= 1:3
%    zii = z(Ind(ii));
%    if(zii>zcenter)
%        if(zii>zMAX2)
%            xHI_matrixTemp(ii,:,:,:) = fftn(ones(N,N,N));
%            %Neut_matrix(ii,:,:,:) =fftn( ones(N,N,N));
%        else
%            load(strcat(pathname_Data2,'xHI_',num2str(zii),'_',num2str(ncube),'_',num2str(fstar),'_',num2str(flag),...
%                '_',num2str(flagM),'_',num2str(XeffTerm),'_',num2str(Ispec),'_',num2str(Reion),'_',num2str(feedback),'_',...
%                num2str(p),'_',num2str(pop),'_',num2str(FSfunc),'_',num2str(photoheatingVersion),'.mat'));% 
%      
%            xHI_matrixTemp(ii,:,:,:) =fftn(max(0,min(1,xHI)));%xe at each pixel at zii
%            xHI = [];
%        end
%    end
%end      
%zextrap = z(Ind);
%xHI_matrixCurrent(:,:,:) =squeeze(interp1(log(1+zextrap(2:3)),xHI_matrixTemp(2:3,:,:,:),log(1+zextrap(1)),'linear','extrap'));
%xHI_matrixTemp=[];
%
%load(strcat(pathname_Data1,'xe_',num2str(zcenter),'_',num2str(ncube),'_',num2str(fstar),'_',num2str(flag),...
%                                '_',num2str(flagM),'_',num2str(XeffTerm),'_',num2str(Ispec),'_',num2str(Reion),'_',num2str(feedback),...
%                                '_',num2str(p),'_',num2str(pop),'_',num2str(FSfunc),'_',num2str(photoheatingVersion),'.mat'));% total LyA flux
%xepix = fftn(xe);%xe at each pixel at zii                           
%xe = [];
%
%
%dfdt_interp = zeros(2,N,N,N);
%dfdt_interpMQ = zeros(2,N,N,N);
%xe_interp = zeros(2,N,N,N);  
%xHI_interp = zeros(2,N,N,N);
%Maxfcoll = zeros(N,N,N);
%eps = zeros(N,N,N);
%Lion = zeros(N,N,N);
%Jalpha=zeros(N,N,N);
%JAX=zeros(N,N,N);
%JA=zeros(N,N,N);
%JLW21=zeros(N,N,N);
%epsMQ = zeros(N,N,N);
%LionMQ = zeros(N,N,N);
%JAXMQ=zeros(N,N,N);
%for ii= 1:Nshells  %integral
%    R=Lpix*(Rmin(ii)+Rmax(ii))/2; % comoving radius of ring
%    zshell = getRtoz(R,zcenter);
%    if(zshell<ceil(zsM)&(zshell>zcenter))
%        izs = find(z>zshell,1,'first');
%        i_interp=0;
%        for iz=izs-1:izs
%            i_interp=i_interp+1;
%            indZ=find(Ind==find(z==z(iz)));
%            dfdt_interp(i_interp,:,:,:) = LWgetGasShell2(squeeze(dfdt_matrix(indZ,:,:,:)),Rmin(ii),Rmax(ii));% df/dt in [sec^-1] units
%            %xe_interp(i_interp,:,:,:) = LWgetGasShell2(squeeze(xe_matrix(indZ,:,:,:)),0,Rmax(ii))/(4*pi*Rmax(ii)^3/3); 
%            %xHI_interp(i_interp,:,:,:) = LWgetGasShell2(squeeze(xHI_matrix(indZ,:,:,:)),0,Rmax(ii))./(4*pi*Rmax(ii)^3/3); 
%            if(MQ)
%                 dfdt_interpMQ(i_interp,:,:,:) = LWgetGasShell2(squeeze(dfdt_matrixMQ(indZ,:,:,:)),Rmin(ii),Rmax(ii));% df/dt in [sec^-1] units
%            end
%         % A added (load history)         
%             if(z(iz)>zMAX2)
%                             xe =  (1.716e-11*z(iz).^3-3.823e-9*z(iz).^2+4.941e-7*z(iz)+0.0001056).*ones(N,N,N);
%                             xe_matrix(:,:,:) =fftn(xe);
%             else
%                             load(strcat(pathname_Data1,'xe_',num2str(z(iz)),'_',num2str(ncube),'_',num2str(fstar),'_',num2str(flag),...
%                            '_',num2str(flagM),'_',num2str(XeffTerm),'_',num2str(Ispec),'_',num2str(Reion),'_',num2str(feedback),...
%                            '_',num2str(p),'_',num2str(pop),'_',num2str(FSfunc),'_',num2str(photoheatingVersion),'.mat'));% total LyA flux
%                            xe_matrix(:,:,:) = fftn(xe);%xe at each pixel at zii
%                            xe = [];
%             end
%                   xe_interp(i_interp,:,:,:) = LWgetGasShell2(squeeze(xe_matrix(:,:,:)),0,Rmax(ii))/(4*pi*Rmax(ii)^3/3); 
%            if(z(iz)>zcenter)
%                        if(z(iz)>zMAX2)
%                            xHI_matrix(:,:,:) = fftn(ones(N,N,N));
%                            %Neut_matrix(ii,:,:,:) =fftn( ones(N,N,N));
%                        else
%                            load(strcat(pathname_Data2,'xHI_',num2str(z(iz)),'_',num2str(ncube),'_',num2str(fstar),'_',num2str(flag),...
%                                '_',num2str(flagM),'_',num2str(XeffTerm),'_',num2str(Ispec),'_',num2str(Reion),'_',num2str(feedback),'_',...
%                                num2str(p),'_',num2str(pop),'_',num2str(FSfunc),'_',num2str(photoheatingVersion),'.mat'));% total LyA flux
%
%                            xHI_matrix(:,:,:) =fftn(max(0,min(1,xHI)));%xe at each pixel at zii
%                            xHI = [];
%                        end
%            end
%              if z(iz)==zcenter
%                        xHI_matrix=xHI_matrixCurrent;
%              end 
%            xHI_interp(i_interp,:,:,:) = LWgetGasShell2(squeeze(xHI_matrix(:,:,:)),0,Rmax(ii))./(4*pi*Rmax(ii)^3/3); 
%            %
%        end
%        xeshell =min(max(xe_grid),max(min(xe_grid),squeeze(interp1(log(1+z(izs-1:izs)),xe_interp,log(1+zshell)))));
%        %mean(mean(mean(xeshell)))
%        xHIshell = min(max(xHI_grid),max(min(xHI_grid),squeeze(interp1(log(1+z(izs-1:izs)),xHI_interp,log(1+zshell)))));
%        %mean(mean(mean(xHIshell)))
%        dfdt =squeeze(exp(interp1(log(1+z(izs-1:izs)),log(abs(dfdt_interp)+1e-26),log(1+zshell))));   %dfdt per shell (including the shell volume)
%        SIM21Testing.saveDebugMat(xeshell,['a_xeshell_',num2str(ii)]);
%        SIM21Testing.saveDebugMat(xHIshell,['a_xHIshell_',num2str(ii)]);
%		SIM21Testing.saveDebugMat(dfdt,['a_dfdt_',num2str(ii)]);
%		SIM21Testing.saveDebugMat(xe_interp,['a_xe_interp_',num2str(ii)]);
%        SIM21Testing.saveDebugMat(xHI_interp,['a_xHI_interp_',num2str(ii)]);
%		SIM21Testing.saveDebugMat(dfdt_interp,['a_dfdt_interp_',num2str(ii)]);
%        JA =  JA+dfdt*LyACoefMat(find(zc==zcenter),ii);% cm-2 sr-1 s-1 Hz-1
%        eps =  eps+(Ispec<6).*XeffTerm*dfdt.*10.^(interp2(log10(xHI_grid+1e-16),log10(xe_grid+1e-16),...
%            log10(squeeze(XCoefMat(find(zc==zcenter),ii,:,:))),log10(squeeze(xHIshell)+1e-16),...
%            log10(squeeze(xeshell)+1e-16)));% eV/sec/baryon 
%        Lion =  Lion+(Ispec<6).*XeffTerm*dfdt.*10.^(interp2(log10(xHI_grid+1e-16),log10(xe_grid+1e-16),...
%            log10(squeeze(LionCoefMat(find(zc==zcenter),ii,:,:))),log10(squeeze(xHIshell)+1e-16),...
%            log10(squeeze(xeshell)+1e-16)));% eV/sec/baryon 
%        
%        JAX = JAX+ (Ispec<6).*XeffTerm*dfdt.*10.^(interp2(log10(xHI_grid+1e-16),log10(xe_grid+1e-16),...
%            log10(squeeze(LyAXCoefMat(find(zc==zcenter),ii,:,:))),log10(squeeze(xHIshell)+1e-16),...
%            log10(squeeze(xeshell)+1e-16)));
%
%        SIM21Testing.saveDebugMat(JA,['a_JA_',num2str(ii)]);
%		SIM21Testing.saveDebugMat(eps,['a_eps_',num2str(ii)]);
%		SIM21Testing.saveDebugMat(Lion,['a_Lion_',num2str(ii)]);
%		SIM21Testing.saveDebugMat(JAX,['a_JAX_',num2str(ii)]);
%
%        %-------------reionization------------------% 
%       if(R<70) 
%            fcoll = LWgetGasShell2(fcoll_matrix,0,Rmax(ii))/(4*pi*Rmax(ii)^3/3);
%            xeR0 = LWgetGasShell2(xepix,0,Rmax(ii))/(4*pi*Rmax(ii)^3/3);       
%            Maxfcoll = max(fcoll+xeR0/zeta,Maxfcoll);
%       end
%       dfdt = [];  
%       fcoll = []; 
%    end
%end
%SIM21Testing.saveDebugMat(JA,'a_JA');
%SIM21Testing.saveDebugMat(eps,'a_eps');
%SIM21Testing.saveDebugMat(Lion,'a_Lion');
%SIM21Testing.saveDebugMat(JAX,'a_JAX');

cd('~/Work/SIM21/ParamStudy/');
%global chosenNode
%chosenNode = 'compute-0-62';
%SIM21Utils.runSimulation(tc1,tc1.name);
%SIM21Utils.runSimulation(tc2,tc2.name);
%SIM21Utils.runSimulation(tc3,tc3.name);
%clear global chosenNode;