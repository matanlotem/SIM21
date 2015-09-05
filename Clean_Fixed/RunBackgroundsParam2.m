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

function RunBackgroundsParam2(outputExt,MyCube,MyStar,MyVBC,MyVc,MyFX,MySED,MyTau,MyFeed,DelayParam,MyPop,FSfunc,photoheatingVersion,varargin)
    addpath('/a/home/cc/tree/taucc/students/physics/matanlotem/Work/SIM21/lib/');
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

    %global pathname
    %pathname = '/scratch300/matanlotem/';

    global vbc_cube
    global delta_cube
    global gridInterpKey;

    global pathname_Data1 % e.g. scratch (to save JLW, x_e, Lion, eps, Jalpha)
    global pathname_Data2 % e. g. scratch300 (to save  xHI, Neut, TK, T21)  
    global pathname_DataBackgrounds % e. g. scratch300 (to save ICs, Window functions, grids)

    if ~ isempty(outputExt)
        if outputExt(end) ~= '/'
            outputExt = [outputExt,'/'];
        end
    end
    pathname_Data1 = ['/scratch/matanlotem/Data/',outputExt];
    pathname_Data2 = ['/scratch300/matanlotem/Data/',outputExt];
    if exist(pathname_Data1) ~= 7
        mkdir(pathname_Data1);
    end
    if exist(pathname_Data2) ~= 7
        mkdir(pathname_Data2);
    end

    pathname_DataBackgrounds = '/scratch300/matanlotem/DataBackgrounds_withPlanck/';
    
    global ID

    photoheatingOn = 0;

    for icube =1:length(MyCube)
        ncube = MyCube(icube)
        delta_main=importdata(strcat(pathname_DataBackgrounds,'my',num2str(ncube),'_d.dat'));
        vbc_main=importdata(strcat(pathname_DataBackgrounds,'my',num2str(ncube),'_v.dat'));
        
        
        delta_cube = delta_main;
        vbc_cube = vbc_main;
        delta = delta_main;
        if photoheatingVersion==1 || photoheatingVersion==2
            photoheatingOn = 1;
        end
        
    %--------------------Simulation------------------------------%
        z = (6:1:120);
        Tau = MyTau;
        pop = MyPop;
        flag = MyVBC;
        flagM = MyVc;
        fstar = MyStar;
        Ispec = MySED;
        XeffTerm = MyFX;
        p = DelayParam;
        feedback = MyFeed;
        %--------preliminary LW (for zeta) -----------------------------
        if (feedback == 1)
            zMAX = 118;
        else
            zMAX = 66;
        end
        
        %---------treshold for reion---------------------------
        Zcube =  0
        delta_Z=importdata(strcat(pathname_DataBackgrounds,'my',num2str(Zcube),'_d.dat'));
        vbc_Z=importdata(strcat(pathname_DataBackgrounds,'my',num2str(Zcube),'_v.dat'));

        delta_cube = delta_Z;  
        vbc_cube = vbc_Z;

        if isempty(cell2mat(varargin))
            zeta=19.48*fstar/0.05;
        else
            zeta=cell2mat(varargin);
        end

        
        %------------LyA, Xrays, T21, Tk, LW ... -----------------------------
        delta_cube = delta_main;  
        vbc_cube = vbc_main;
        
        ID = SIM21Utils.getID(MyCube,MyStar,MyVBC,MyVc,MyFX,MySED,MyTau,MyFeed,DelayParam,MyPop,FSfunc,photoheatingVersion);
        gridInterpKey = [ID,'.mat'];
        
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
    vbc_cube = [];
    delta_cube = [];
    JLW21 = [];


    % p = profile('info');
    % profile summary report  
    toc;
end
