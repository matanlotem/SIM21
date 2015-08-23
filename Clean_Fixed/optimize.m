pathname_Data1 = '/scratch/matanlotem/Data/';
pathname_Data2 = '/scratch300/matanlotem/Data/';
pathname_DataBackgrounds = '/scratch300/matanlotem/DataBackgrounds_withPlanck/';

zcenter = 60;
ncube = 0;
fstar = 0.05;
flag = 1;
flagM = 16.5;
XeffTerm = 1;
Ispec = 1;
Tau = 0.075;
zeta = 20;
feedback = 0;
p = 0;
pop = 2;
FSfunc = 1;
photoheatingOn = 1;
photoheatingVersion = 2;

ID = '_0_0.05_1_16.5_1_1_0.075_0_0_2_1_2';

delta_cube = importdata(strcat(pathname_DataBackgrounds,'my',num2str(ncube),'_d.dat'));

OBackgroundsParamII(zcenter,ncube,fstar,flag,flagM,XeffTerm,Ispec,Tau,zeta,feedback,p,pop,FSfunc,photoheatingOn,photoheatingVersion);