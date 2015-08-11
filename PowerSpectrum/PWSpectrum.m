function ID = ParamStudy(MyCube,MyStar,MyVBC,MyVc,MyFX,MySED,MyTau,MyFeed,DelayParam,MyPop,FSfunc,photoheatingVersion,zeta) 
tic;
global pathname_Data1 % e.g. scratch (to save JLW, x_e, Lion, eps, Jalpha)
pathname_Data1 = '/scratch/matanlotem/Data/';

ID = ['_' num2str(MyCube)...
      '_' num2str(MyStar) '_' num2str(MyVBC) '_' num2str(MyVc)...
      '_' num2str(MyFX) '_' num2str(MySED) '_' num2str(MyTau)...
      '_' num2str(MyFeed) '_' num2str(DelayParam) '_' num2str(MyPop) '_' num2str(FSfunc) '_' num2str(photoheatingVersion)];


global vbc_cube
global delta_cube
delta_cube=importdata(strcat('/scratch300/matanlotem/DataBackgrounds_withPlanck/my',num2str(MyCube),'_d.dat'));
vbc_cube=importdata(strcat('/scratch300/matanlotem/DataBackgrounds_withPlanck/my',num2str(MyCube),'_v.dat'));
N=length(vbc_cube);
Lpix=3;
ep=10;
Lx=Lpix*N;
Nmu = 11;

Nvec = MyCube;
%zspec = [6:0.1:15, 16:1:40];
zspec = [6:0.1:6.2];
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
            [Pk,Pk_iso,Pk_X,Pk_del,Pmu4,Pmu2,Pmu0,K,T21] = getPKIIRec(zspec(indz),MyCube,MyStar,MyVBC,MyVc,MyFX,MySED,MyTau,MyFeed,DelayParam,MyPop,FSfunc,photoheatingVersion); 
            PowerMat(indz,:)= Pk;
            PowerMat_iso(indz,:)= Pk_iso;
            PowerMat_X(indz,:)= Pk_X; 
            PowerMat_del(indz,:)= Pk_del;
            PowerMat_4(indz,:)= Pmu4;
            PowerMat_2(indz,:)= Pmu2;
            PowerMat_0(indz,:)= Pmu0;
            
            [Tb,T21] = getTbcube(zspec(indz),MyCube,MyStar,MyVBC,MyVc,MyFX,MySED,MyTau,MyFeed,DelayParam,MyPop,FSfunc,photoheatingVersion);
            del = real(max(-0.9,min(1,delta_cube*LWgetDz(zspec(indz))/LWgetDz(40))));
            delta_v = getDV(del,Lx);
            Tlin = Tb.*(1+delta_v); 
            T21 = mean(mean(mean(Tlin)));
            Tlin = real((Tlin-T21)./T21);
            [Pk,Pkm,Kout,Muout,nk,nkm] = getPkMu(Tlin,Lx,ep,Nmu);
            Pkm = Pkm*T21.^2;
            PowerMatMu(indz,:,:) = Pkm;
end

save(strcat('/scratch300/matanlotem/ParamStudy/TPowerMat',ID,'.mat'),'PowerMat');
save(strcat('/scratch300/matanlotem/ParamStudy/TPowerMat_iso',ID,'.mat'),'PowerMat_iso');
save(strcat('/scratch300/matanlotem/ParamStudy/TPowerMat_X',ID,'.mat'),'PowerMat_X');
save(strcat('/scratch300/matanlotem/ParamStudy/TPowerMat_del',ID,'.mat'),'PowerMat_del');
save(strcat('/scratch300/matanlotem/ParamStudy/TPowerMat_0',ID,'.mat'),'PowerMat_0');
save(strcat('/scratch300/matanlotem/ParamStudy/TPowerMat_2',ID,'.mat'),'PowerMat_2');
save(strcat('/scratch300/matanlotem/ParamStudy/TPowerMat_4',ID,'.mat'),'PowerMat_4');
save(strcat('/scratch300/matanlotem/ParamStudy/TPowerMat_mu',ID,'.mat'),'PowerMatMu');
toc;
end