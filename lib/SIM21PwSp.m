classdef SIM21PwSp
    methods(Static)
        %%%function ID = getPwSp(MyCube,MyStar,MyVBC,MyVc,MyFX,MySED,MyTau,MyFeed,DelayParam,MyPop,FSfunc,photoheatingVersion,zeta) 
        function ID = getPwSp(runCase) 
            tic;
            global pathname_Data1 % e.g. scratch (to save JLW, x_e, Lion, eps, Jalpha)
            global pathname_Data2
            global ID
            global vbc_cube
            global delta_cube

            pathname_Data1 = runCase.tmpDataPath;
            pathname_Data2 = runCase.dataPath;
            pathname_DataBackgrounds = SIM21Utils.paths.dataBackgrounds;
            pathname_Output = runCase.outputPath;
            ID = runCase.ID;

            delta_cube = importdata([pathname_DataBackgrounds,'my',num2str(runCase.ncube),'_d.dat']);
            vbc_cube = importdata([pathname_DataBackgrounds,'my',num2str(runCase.ncube),'_v.dat']);
            
            N=length(vbc_cube);
            Lpix=3;
            ep=10;
            Lx=Lpix*N;
            Nmu = 11;

            Nvec = runCase.ncube;
            zspec = [6:0.1:15, 16:1:40];
            %zspec = 6:0.1:6.2; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            PowerMat=zeros(length(zspec),100);
            PowerMat_iso=zeros(length(zspec),100);
            PowerMat_X=zeros(length(zspec),100);
            PowerMat_del=zeros(length(zspec),100);
            PowerMat_0=zeros(length(zspec),100);
            PowerMat_2=zeros(length(zspec),100);
            PowerMat_4=zeros(length(zspec),100);
            PowerMatMu=zeros(length(zspec),Nmu,100);
            
            for indz = 1:length(zspec)
                disp(['z = ',num2str(zspec(indz))]);
                PS=[];  
                Pk_iso=[];
                Pk_X=[];
                Pk_del=[];
                %[Pk,Pk_iso,Pk_X,Pk_del,Pmu4,Pmu2,Pmu0,K,T21] = getPKIIRec(zspec(indz),MyCube,MyStar,MyVBC,MyVc,MyFX,MySED,MyTau,MyFeed,DelayParam,MyPop,FSfunc,photoheatingVersion); 
                %%%%%
                
                del = real(max(-0.9,min(1,delta_cube*LWgetDz(zspec(indz))/LWgetDz(40))));
                %%%[Tb,T21] = SIM21PwSp.getTbcube(zspec(indz),MyCube,MyStar,MyVBC,MyVc,MyFX,MySED,MyTau,MyFeed,DelayParam,MyPop,FSfunc,photoheatingVersion);
                disp('getTbcube');
                [Tb,T21] = SIM21PwSp.getTbcube(zspec(indz),runCase);
                Tb2 = real((Tb-T21)./T21);
                
                disp('getPkRec');
                [Pk,Pk_iso,Pk_X,Pk_del,Pmu4,Pmu2,Pmu0,K,nk] = SIM21PwSp.getPkRec(Tb2,del,Lx,ep);
                Pk = Pk.*T21^2; 
                Pk_iso = Pk_iso.*T21^2;  
                Pk_X = Pk_X.*T21^2;   
                Pk_del = Pk_del.*T21^2;   
                Pmu4 =  Pmu4.*T21^2;
                Pmu2 =  Pmu2.*T21^2;
                Pmu0 =  Pmu0.*T21^2;
                
                %%%%%
                
                PowerMat(indz,:)= Pk;
                PowerMat_iso(indz,:)= Pk_iso;
                PowerMat_X(indz,:)= Pk_X; 
                PowerMat_del(indz,:)= Pk_del;
                PowerMat_4(indz,:)= Pmu4;
                PowerMat_2(indz,:)= Pmu2;
                PowerMat_0(indz,:)= Pmu0;
                
                
                %[Tb,T21] = getTbcube(zspec(indz),MyCube,MyStar,MyVBC,MyVc,MyFX,MySED,MyTau,MyFeed,DelayParam,MyPop,FSfunc,photoheatingVersion);
                %del = real(max(-0.9,min(1,delta_cube*LWgetDz(zspec(indz))/LWgetDz(40))));
                
                %%%% Don't get Pkm
                %%%% disp('getDV');
                %%%% delta_v = SIM21PwSp.getDV(del,Lx);
                %%%% Tlin = Tb.*(1+delta_v); 
                %%%% T21 = mean(mean(mean(Tlin)));
                %%%% Tlin = real((Tlin-T21)./T21);
                %%%% 
                %%%% %[Pk,Pkm,Kout,Muout,nk,nkm] = SIM21PwSp.getPkMu(Tlin,Lx,ep,Nmu);
                %%%% disp('gePktMu');
                %%%% [Pkm,Kout,Muout,nkm] = SIM21PwSp.getPkMu(Tlin,Lx,ep,Nmu);
                %%%% Pkm = Pkm*T21.^2;
                %%%% PowerMatMu(indz,:,:) = Pkm;
            end

            save([pathname_Output,'PowerMat',ID,'.mat'],'PowerMat');
            save([pathname_Output,'PowerMat_iso',ID,'.mat'],'PowerMat_iso');
            save([pathname_Output,'PowerMat_X',ID,'.mat'],'PowerMat_X');
            save([pathname_Output,'PowerMat_del',ID,'.mat'],'PowerMat_del');
            save([pathname_Output,'PowerMat_0',ID,'.mat'],'PowerMat_0');
            save([pathname_Output,'PowerMat_4',ID,'.mat'],'PowerMat_4');
            save([pathname_Output,'PowerMat_2',ID,'.mat'],'PowerMat_2');
            %%%%save([pathname_Output,'N_TPowerMat_mu',ID,'.mat'],'PowerMatMu');
            toc;
        end
        
        
        function [kk,Vpix,mu] = getPkVars(cube,Lx)
            N = length(cube);
            dx = Lx/N;
            pix = dx;
            Vpix = dx^3;
            
            kx = zeros(size(cube));
            ky = zeros(size(cube));
            kz = zeros(size(cube));
            
            for Ind = 0:N-1
                kx(Ind+1,:,:) = -(Ind>(N-1)/2)*(2*pi*(N-Ind))/(N*pix)*ones(1,N,N)+(Ind<=(N-1)/2)*(2*pi*Ind/(N*pix))*ones(1,N,N);
                ky(:,Ind+1,:) = -(Ind>(N-1)/2)*(2*pi*(N-Ind))/(N*pix)*ones(N,1,N)+(Ind<=(N-1)/2)*(2*pi*Ind/(N*pix))*ones(N,1,N); 
                kz(:,:,Ind+1) = -(Ind>(N-1)/2)*(2*pi*(N-Ind))/(N*pix)*ones(N,N,1)+(Ind<=(N-1)/2)*(2*pi*Ind/(N*pix))*ones(N,N,1);
            end
            
            kk = (kx.^2+ky.^2+kz.^2).^0.5;
            
            %--- adding velocity term (20.01.2014) --------
            mu = kz./kk;
            %Ind = find(kk==0); %Ind = kk==0; %%% I think this is a bug - in getPkRec line 30
            Ind = kk==0;
            mu(Ind) = 0;
        end
        
        
        function [KN,Kout] = getKout(Lx,ep)
            k0=2*pi*((2*ep+1)/(2*ep))/Lx;
            powK=0:100;
            KN = k0*((2*ep+1)/(2*ep-1)).^powK;
            
            Kout = (2*pi/Lx)*((2*ep+1)/(2*ep-1)).^(1:100);
        end
        
        
        function [delta_out] = getDV(del,Lx)
            [kk,Vpix,mu] = SIM21PwSp.getPkVars(del,Lx);
            kdel = fftn(del);
            kcube_del = kdel.*mu.^2;
            delta_out = (ifftn(kcube_del));
        end
        
        
        %%%function [Tb,T21] = getTbcube(zii,ncube,fstar,flag,flagM,XeffTerm,Ispec,Reion,feedback,p,pop,FSfunc,photoheatingVersion)
        function [Tb,T21] = getTbcube(zii,runCase)
            global pathname_Data2
            global ID
            global delta_cube
            load(SIM21Utils.getMatrixPath('Planck_parameters'));
            
            del = real(max(-0.9,min(1,delta_cube*LWgetDz(zii)/LWgetDz(40))));
            if(runCase.ncube<100)
                Lpix = 3;
            else
                Lpix = 6;
            end
            
            N=length(delta_cube);
            ep=10;
            Lx=Lpix*N;
            C21 = 27*sqrt(0.15/Om/h^2)*(Ob*h^2/0.023);
            z21cm=5:50;
            zion = [5:0.1:15, 16:1:50];
            
            Ind1 = find(z21cm>zii,1,'first');  
            zint = z21cm(Ind1-1:Ind1);  
            F_interp = zeros(2,N,N,N);
            
            for iz=1:2
                disp(['load ',num2str(iz)]);
                load(strcat(pathname_Data2,'TK_',num2str(zint(iz)),ID,'.mat'));
                disp('getTs');
                Ts = SIM21Gets.getTs(TK,zint(iz),runCase.ncube,runCase.fstar,runCase.vbc,runCase.vc,runCase.fx,runCase.sed,runCase.zeta,runCase.feedback,runCase.delayParam,runCase.pop,runCase.fsfunc,runCase.phVersion);
                %Ts = TK;
                disp('gotTs');
                Ts(Ts==0)=1e-20;   
                deltaTerm = (ones(size(delta_cube))+max(-0.9,min(1,delta_cube*LWgetDz(zint(iz))/LWgetDz(40))));
                
                Tb = C21*sqrt((1+zint(iz))/10).*deltaTerm.*(Ts-2.725*(1+zint(iz)))./Ts;
                F_interp(iz,:,:,:) = Tb;
                Tb =[];        
            end
            
            Tb = squeeze((interp1(log(1+zint),F_interp,log(1+zii)))); 
            F_interp=[];
            
            load([pathname_Data2,'xHI_',num2str(zii),ID,'.mat']);
            load([pathname_Data2,'Neut_',num2str(zii),ID,'.mat']);
            
            Tb =  Tb.*xHI.*Neut;
            T21 = mean(mean(mean(Tb)));
           
        end
        
        function [Pk,Pk_iso,Pk_X,Pk_del,Pmu4,Pmu2,Pmu0,Kout,nk] = getPkRec(cube,del,Lx,ep)
            [kk,Vpix,mu] = SIM21PwSp.getPkVars(cube,Lx);
            [KN,Kout] = SIM21PwSp.getKout(Lx,ep);
            
            kdel = fftn(del);
            kcube = fftn(cube);
            kcube_del = kdel.*mu.^2;
            kcube_new = kcube+kcube_del;
            
            [NNN, bink] = histc(kk,KN);
            
            lKN1=length(KN)-1;
            NMbink = repmat(bink,[1,1,1,lKN1]) == permute(repmat((1:lKN1)',[1,size(bink)]),[2,3,4,1]);
            nk = squeeze(sum(sum(sum(NMbink))))';
            
            function out = calcThis(tmp_cube)
                M = repmat(tmp_cube,[1,1,1,lKN1]).*NMbink;
                out = 1./(Lx*Lx*Lx)*(Vpix^2)*squeeze(sum(sum(sum(M))))'./nk;
                out(nk==0) = 0;
            end
            
            Pk = calcThis(kcube_new.*conj(kcube_new));
            Pk_iso = calcThis(kcube.*conj(kcube));
            Pk_X = calcThis((kcube.*conj(kcube_del))+(kcube_del.*conj(kcube)));
            Pk_del = calcThis(kcube_del.*conj(kcube_del));
            
            Pmu4 = calcThis(kdel.*conj(kdel));
            Pmu2 = calcThis((kcube.*conj(kdel)) + (kdel.*conj(kcube)));
            Pmu0 = Pk_iso;
        end
        
        
        %function [Pk,Pkm,Kout,Muout,nk,nkm] = getPkMu(cube,Lx,ep,Nmu)
        function [Pkm,Kout,Muout,nkm] = getPkMu(cube,Lx,ep,Nmu)
            [kk,Vpix,mu] = SIM21PwSp.getPkVars(cube,Lx);
            [KN,Kout] = SIM21PwSp.getKout(Lx,ep);
            
            Muvec = linspace(-1,1,Nmu+1);
            Muout = (Muvec(1:end-1)+Muvec(2:end))/2;
            
            kcube = fftn(cube);
            
            [NNN, bink] = histc(kk,KN);
            [MMM, mink] = histc(mu,Muvec);
            
            lKN1=length(KN)-1;
            NMbink = repmat(bink,[1,1,1,lKN1]) == permute(repmat((1:lKN1)',[1,size(bink)]),[2,3,4,1]);
            NMmink = repmat(mink,[1,1,1,Nmu]) == permute(repmat((1:Nmu)',[1,size(mink)]),[2,3,4,1]);
            
            %%%%%!!!!!!!!!!!! REQUIRES A LOT OF MEMORY !!!!!!!!!!!!!
            NMbm = repmat(NMbink,[1,1,1,1,Nmu]).*permute(repmat(NMmink,[1,1,1,1,lKN1]),[1,2,3,5,4]);
            
            %%%%%% WHY? This is calculated elsewhere...
            %nk = squeeze(sum(sum(sum(NMbink))))';
            %Pk = 1./(Lx*Lx*Lx)*(Vpix^2)*squeeze(sum(sum(sum(repmat((kcube.*conj(kcube)),[1,1,1,lKN1]).*NMbink))))'./nk;
            %Pk(nk==0) = 0;
            %%%%%%
            NMbink = [];
            NMmink = [];
            
            nkm = squeeze(sum(sum(sum((NMbm)))))';
            NMbmk = squeeze(sum(sum(sum(repmat(kcube.*conj(kcube),[1,1,1,lKN1,Nmu]).*NMbm))))';
            
            Pkm = 1./(Lx*Lx*Lx)*(Vpix^2)*NMbmk./nkm;
            Pkm(nkm==0) = 0;
        end
        
        
% -------------------- Testing -------------------- %        
% ------------------------------------------------- %
        function compareOutputs()
            pathname_Output = '/scratch300/matanlotem/ParamStudy/';
            TFiles = dir([pathname_Output,'P*']);
            for i = 1:length(TFiles)
                TFile = TFiles(i);
                flag = isequal(importdata([pathname_Output,TFile.name]),importdata([pathname_Output,'N_',TFile.name]));
                disp([num2str(flag),' ',TFile.name]);
            end
        end
        
        
        function testParams()
            SIM21PwSp.getPwSp(0,0.05,1,16.5,1,1,0.075,0,0,2,1,2,20);
        end
    end
end