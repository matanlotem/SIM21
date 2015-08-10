% Run BackgroundsReionII.m runs function BackgroundsReionII.m for each redshift 

% Input parameters: 
% 1. MyCube: 0 to 19 selects a set of initial conditions. Can be a vector
% 2. MyStar: star formation efficiency (e.g. = 0.1). Separate efficiencies
% for molecular cooling and atomic cooling (should be a vector of the same
% length as MyVc, specify fstar for every Vcircular). 
% 3. MyVBC: 1/0, with/without vbc 
% 4. MyVc: the value of circular velocity (4.2 or 3.7 km/s for molecular cooling, 16.5 km/s atomic cooling, 35 km/s heavy halos)
% 5. MyFX: Normalization of Xrays (f_X), XeffTerm = 1 corresponds to 3e40
% [erg s/yr/M_sun]. Can be a vector.
% 6. MySED:  1 hard X-ray spectrum. 0 is for the 50% mixture (geometric mean), 2 is for the old case of power-law. 
% Can be a vector, 3 50% and Mini Qasars, 4 new spectrum and MQ, 5 old  %Q
% spectrum and MQ, 6 only MQ.       %Q
% 7. MyTau: the value of total optical depth of CMB (tau = 0.089 - early
% reionization, tau = 0.075 - late reionization). Can be a vector
% 8. MyFeed: 0/1 - without/with LW, can be a vector
% 9. DelayParam: any from 0 to 1 (e.g. 0.75/0.5 - strong/weak feedback).
% Vector of the same length as MyFeed
% 10. MyPop: 2/3 determines stellar population. For PopII use 
% f_starII =  3*f_starIII

% Cubes of LyA, Xrays, LW, ionizing radiative backgrounds, T21cm, Tgas, xe
% and neutral fraction at each redshift are saved to memory

function RunBackgroundsParam(MyCube,MyStar,MyVBC,MyVc,MyFX,MySED,MyTau,MyFeed,DelayParam,MyPop,FSfunc,photoheatingVersion,varargin) 
tic;
% MyCube = 0;
% MyStar = 0.05;
% MyVBC=1;
% MyVc =4.2;
% MyFX=1;
% MySED=1;
% MyTau=0.075;
% MyFeed=1;
% DelayParam=0.5;
% MyPop = 3;
 
   % profile on -history;   

global pathname
pathname = '/scratch300/matanlotem/';

global vbc_cube
global delta_cube
global gridInterpKey;

global pathname_Data1 % e.g. scratch (to save JLW, x_e, Lion, eps, Jalpha)
pathname_Data1 = '/scratch/matanlotem/Data/';
global pathname_Data2 % e. g. scratch300 (to save  xHI, Neut, TK, T21)  
pathname_Data2 = '/scratch300/matanlotem/Data/';
global pathname_Data3 % e. g. scratch300 (to save ICs, Window functions, grids)  
pathname_Data3 = '/scratch/matanlotem/';

photoheatingOn = 0;

for icube =1:length(MyCube)
   ncube = MyCube(icube)
    delta_main=importdata(strcat(pathname,'DataBackgrounds_withPlanck/my',num2str(ncube),'_d.dat'));
    vbc_main=importdata(strcat(pathname,'DataBackgrounds_withPlanck/my',num2str(ncube),'_v.dat'));
    
    
    delta_cube = delta_main;
    vbc_cube = vbc_main;
	delta = delta_main;
    if photoheatingVersion==1 || photoheatingVersion==2
                                photoheatingOn = 1;
                                end
    
%--------------------Simulation------------------------------%
    z = (6:1:120);
    for indt = 1:length(MyTau)
    	Tau = MyTau(indt)
        for indp = 1:length(MyPop)
            pop = MyPop(indp)
        for indv = 1:length(MyVBC)
                flag = MyVBC(indv);
                for indM = 1:length(MyVc)
                    flagM = MyVc(indM);
                    fstar = MyStar(indM);
                    i=0;
                  % zeta = getZeta2(Tau,flag,flagM) % the value of zeta for reionization 
                 
                   for indsed = 1:length(MySED)
                    	Ispec = MySED(indsed)
                        for indfx = 1:length(MyFX)
                        	XeffTerm = MyFX(indfx)
                            for indf = 1:length(MyFeed)
                                p = DelayParam(indf);
                                feedback = MyFeed(indf); 
                                %--------preliminary LW (for zeta) -----------------------------
                                if (feedback == 1)
                                    zMAX = 118;
                                    i=0;
%                                     while(z(find(z==zMAX)-i+1)>6)
%                                          zcenter = z(find(z==zMAX)-i)
%                                          getLW(zcenter,ncube,fstar,flag,flagM,XeffTerm,Ispec,Tau,feedback,p,pop,FSfunc,photoheatingOn,photoheatingVersion);   
%                                          i=i+1;
%                                     end
                                else
                                    zMAX = 66;
                                end
                                
                                %---------treshold for reion---------------------------
                                Zcube =  0
                                delta_Z=importdata(strcat(pathname,'DataBackgrounds_withPlanck/my',num2str(Zcube),'_d.dat'));
                                vbc_Z=importdata(strcat(pathname,'DataBackgrounds_withPlanck/my',num2str(Zcube),'_v.dat'));
    
                                %delta_Z=importdata(strcat('/scratch/prun/P_DATA/my',num2str(Zcube),'_d.dat'));
                                %vbc_Z=importdata(strcat('/scratch/prun/P_DATA/my',num2str(Zcube),'_v.dat'));
                                delta_cube =delta_Z;  
                                vbc_cube = vbc_Z;
                                    %zeta = getZeta2(Zcube,Tau,flag,flagM,fstar,feedback,p,pop)
                                if isempty(cell2mat(varargin))
                                    zeta=19.48*fstar/0.05;
                                else
                                    zeta=cell2mat(varargin);
                                end

                                
                                %------------LyA, Xrays, T21, Tk, LW ... -----------------------------
                                delta_cube =delta_main;  
                                vbc_cube = vbc_main;
                                
                                gridInterpKey = ['_' num2str(ncube)...
                                     '_' num2str(fstar) '_' num2str(flag) '_' num2str(flagM)...
                                     '_' num2str(XeffTerm) '_' num2str(Ispec) '_' num2str(Tau)...
                                     '_' num2str(feedback) '_' num2str(p) '_' num2str(pop) '_' num2str(FSfunc) '_' num2str(photoheatingVersion) '.mat'];
                                
                                                               
                                i=0;
                                if photoheatingOn==0
                                 while(z(find(z==zMAX)-i+1)>6)
                                     zcenter = z(find(z==zMAX)-i)
                                     BackgroundsParamII(zcenter,ncube,fstar,flag,flagM,XeffTerm,Ispec,Tau,zeta,feedback,p,pop,FSfunc,photoheatingOn,photoheatingVersion); 
                                     %BackgroundsPix(zcenter,ncube,fstar,flag,flagM,XeffTerm,Ispec,Tau,zeta,feedback,p);         
                                     i=i+1;
                                 end
                                %high resolution for reionization 
                                zION = 6:0.1:15;
                                for indz=1:length(zION)
                                    zi = zION(indz)
                                    xHImean = BackgroundsParamIIRes(zi,ncube,fstar,flag,flagM,XeffTerm,Ispec,Tau,zeta,feedback,p,pop,FSfunc,photoheatingOn,photoheatingVersion); 
                                    %xHImean = BackgroundsPixRes(zi,ncube,fstar,flag,flagM,XeffTerm,Ispec,Tau,zeta,feedback,p);
                                end
                                %-------------------------------------
                                end
                                
                                if photoheatingOn==1
                                while(z(find(z==zMAX)-i+1)>6)
                                     zcenter = z(find(z==zMAX)-i)
                                     BackgroundsParamII(zcenter,ncube,fstar,flag,flagM,XeffTerm,Ispec,Tau,zeta,feedback,p,pop,FSfunc,photoheatingOn,photoheatingVersion); 
                                     %BackgroundsPix(zcenter,ncube,fstar,flag,flagM,XeffTerm,Ispec,Tau,zeta,feedback,p);         
                                        if zcenter<16 && zcenter>6
                                        for k=1:9
                                         zres=zcenter-0.1*k
                                         xHImean = BackgroundsParamIIResPH(zres,ncube,fstar,flag,flagM,XeffTerm,Ispec,Tau,zeta,feedback,p,pop,FSfunc,photoheatingOn,photoheatingVersion); 
                                         clear grid_interpSF2;
                                        end
                                        end
                                     i=i+1;
                                end
                                end
                                 
                            end
                        end
                    end
                end
        end
    end
    end
end
vbc_cube = [];
delta_cube = [];
JLW21 = [];



% p = profile('info');
% profile summary report  
toc;
end
