addpath('/a/home/cc/tree/taucc/students/physics/matanlotem/Work/SIM21/lib/');

c = cases55(3);
z = 9;


global pathname
global pathname_Data1
global pathname_Data2
global pathname_DataBackgrounds
global delta_cube

pathname = '/scratch300/matanlotem/';
pathname_Data1 = c.tmpDataPath;
pathname_Data2 = c.dataPath;
pathname_DataBackgrounds = SIM21Utils.paths.dataBackgrounds;
delta_cube=importdata(strcat(pathname_DataBackgrounds,'my',num2str(c.ncube),'_v.dat'));

zcenter = z;

N = 128;
zMAX2 = 60;
ID = c.ID;
load([pathname_Data1,'Lion_',num2str(zcenter),ID,'.mat']);
load([pathname_Data1,'eps_',num2str(zcenter),ID,'.mat']);
load([pathname_Data1,'xe_',num2str(zcenter),ID,'.mat']);
load([pathname_Data2,'TK_',num2str(zcenter),ID,'.mat']);
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
A1 = A0+SIM21Gets.getf(zcenter,A0,B0,Lion0).*(x1-x0);
B1 = B0+SIM21Gets.getg(zcenter,A0,B0,Lion0,eps0).*(x1-x0);

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
        load([pathname_Data1,'Lion_',num2str(zi(indz)),ID,'.mat'],'Lion');
        load([pathname_Data1,'eps_',num2str(zi(indz)),ID,'.mat'],'eps');

        Lion_extrap(indz,:,:,:) = Lion;% df/dt in [sec^-1] units 
        eps_extrap(indz,:,:,:) = eps;% df/dt in [sec^-1] units    
    end
end
Lion1 = squeeze(interp1(log10(1+zi), Lion_extrap,log10(1+z1),'linear','extrap'));% extrap Lion at the intermediate step
eps1 =  squeeze(interp1(log10(1+zi), eps_extrap,log10(1+z1),'linear','extrap'));% extrap Lion at the intermediate step
Lion_extrap = [];
eps_extrap=[];
%-----step 3

A2 = A0+ SIM21Gets.getf(z1,A1,B1,Lion1).*(x2-x0);
B2 = B0+ SIM21Gets.getg(z1,A1,B1,Lion1,eps1).*(x2-x0);

xe = exp(A2)-1;
TK = exp(B2);

%save([pathname_Data1,'xe_',num2str(zcenter-1),ID,'.mat'],'xe');
%save([pathname_Data2,'TK_',num2str(zcenter-1),ID,'.mat'],'TK');
%

%  count = 4                  
%TK = [];
%xe = [];
%Lion1=[];
%eps1=[];
%Lion = [];
%eps = [];

%maxs = cat(1,p.workCases,zeros(1,length(p.workCases)));
%for i=1:length(maxs)
%    caseNum = maxs(1,i)
%    maxs(2,i)=max(max(max(max(FdebugNaNs.loadMat(SIM21Analysis.XeMagic,SIM21Analysis.XeZ,p.tmpDataPath,p.paramCases(caseNum))))));
%end