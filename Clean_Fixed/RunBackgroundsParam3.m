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

function RunBackgroundsParam2(dataPath,tmpDataPath,MyCube,MyStar,MyVBC,MyVc,MyFX,MySED,MyTau,MyFeed,DelayParam,MyPop,FSfunc,photoheatingVersion,zeta,varargin)
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
    global ID

    pathname_DataBackgrounds = '/scratch300/matanlotem/DataBackgrounds_withPlanck/';
    pathname_Data1 = tmpDataPath;
    pathname_Data2 = dataPath;
    if exist(pathname_Data1) ~= 7
        mkdir(pathname_Data1);
    end
    if exist(pathname_Data2) ~= 7
        mkdir(pathname_Data2);
    end

    ncube = MyCube;
    Tau = MyTau;
    pop = MyPop;
    flag = MyVBC;
    flagM = MyVc;
    fstar = MyStar;
    Ispec = MySED;
    XeffTerm = MyFX;
    p = DelayParam;
    feedback = MyFeed;

    
    delta_cube=importdata(strcat(pathname_DataBackgrounds,'my',num2str(ncube),'_d.dat'));
    vbc_cube=importdata(strcat(pathname_DataBackgrounds,'my',num2str(ncube),'_v.dat'));

    ID = SIM21Utils.getID(ncube,fstar,flag,flagM,XeffTerm,Ispec,Tau,feedback,p,pop,FSfunc,photoheatingVersion)
    zeta
    gridInterpKey = [ID,'.mat'];

    photoheatingOn = 0;
    if photoheatingVersion==1 || photoheatingVersion==2
        photoheatingOn = 1;
    end
    
%--------------------Simulation------------------------------%
    %z = (6:1:120);
    z = [120:-1:6];
    %--------preliminary LW (for zeta) -----------------------------
    if (feedback == 1)
        zMAX = 118;
    else
        zMAX = 61;
    end
    if ~isempty(cell2mat(varargin))
        zMAX = min(varargin{1},zMAX);
    end
    z = z(find(z==zMAX):end);
    

    %------------LyA, Xrays, T21, Tk, LW ... -----------------------------
    if photoheatingOn==0
        for zcenter = z
            disp(zcenter);
            BackgroundsParamII(zcenter,ncube,fstar,flag,flagM,XeffTerm,Ispec,Tau,zeta,feedback,p,pop,FSfunc,photoheatingOn,photoheatingVersion); 
        end
        %high resolution for reionization 
        zION = 6:0.1:15;
        for zi = zION
            disp(zION);
            xHImean = BackgroundsParamIIRes(zi,ncube,fstar,flag,flagM,XeffTerm,Ispec,Tau,zeta,feedback,p,pop,FSfunc,photoheatingOn,photoheatingVersion); 
        end
        %-------------------------------------
    end

    if photoheatingOn==1
        for zcenter = z
            disp(zcenter);
            BackgroundsParamII(zcenter,ncube,fstar,flag,flagM,XeffTerm,Ispec,Tau,zeta,feedback,p,pop,FSfunc,photoheatingOn,photoheatingVersion); 
            if zcenter<16 && zcenter>6
                for zres = [zcenter-0.1:-0.1:zcenter-0.9]
                    disp(zres);
                    xHImean = BackgroundsParamIIResPH(zres,ncube,fstar,flag,flagM,XeffTerm,Ispec,Tau,zeta,feedback,p,pop,FSfunc,photoheatingOn,photoheatingVersion); 
                    clear grid_interpSF2;
                end
            end
        end
    end

    vbc_cube = [];
    delta_cube = [];
    JLW21 = [];

    toc;
end
