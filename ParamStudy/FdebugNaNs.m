classdef FdebugNaNs
    methods(Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Running Stuff Tools
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function runSimulation_Old(c)
            disp('==Running Full Simulation==');
            
            addpath('/a/home/cc/tree/taucc/students/physics/matanlotem/Work/SIM21/Clean/');
            RunBackgroundsParam(c.ncube,c.fstar,c.vbc,c.vc,c.fx,c.sed,c.tau,c.feedback,c.delayParam,c.pop,c.fsfunc,c.photoHeating,c.zeta);
            rmpath('/a/home/cc/tree/taucc/students/physics/matanlotem/Work/SIM21/Clean/');
        end
        
        
        function runZ_Old(z,c)
            global pathname
            global pathname_Data1
            global pathname_Data2
            global delta_cube
            
            disp(['==Running z=',num2str(z),'==']);
            tic;
            
            addpath('/a/home/cc/tree/taucc/students/physics/matanlotem/Work/SIM21/Clean/');
            delta_cube=importdata(strcat(pathname,'DataBackgrounds_withPlanck/my',num2str(c.ncube),'_d.dat'));
            BackgroundsParamII(z,c.ncube,c.fstar,c.vbc,c.vc,c.fx,c.sed,c.tau,c.zeta,c.feedback,c.delayParam,c.pop,c.fsfunc,~~c.photoHeating,c.photoHeating);
            rmpath('/a/home/cc/tree/taucc/students/physics/matanlotem/Work/SIM21/Clean/');
            
            toc;
        end
        
        
        function runSimulation(c)
            disp('==Running Full Simulation==');
            
            addpath('/a/home/cc/tree/taucc/students/physics/matanlotem/Work/SIM21/Clean_Fixed/');
            RunBackgroundsParam(c.ncube,c.fstar,c.vbc,c.vc,c.fx,c.sed,c.tau,c.feedback,c.delayParam,c.pop,c.fsfunc,c.photoHeating,c.zeta);
            rmpath('/a/home/cc/tree/taucc/students/physics/matanlotem/Work/SIM21/Clean_Fixed/');
        end
        
        
        function runZ(z,c)
            global pathname
            global pathname_Data1
            global pathname_Data2
            global pathname_DataBackgrounds
            global delta_cube
            global ID
            
            disp(['==Running z=',num2str(z),'==']);
            tic;
            
            addpath('/a/home/cc/tree/taucc/students/physics/matanlotem/Work/SIM21/Clean_Fixed/');
            ID = c.ID;
            delta_cube=importdata(strcat(pathname,'DataBackgrounds_withPlanck/my',num2str(c.ncube),'_d.dat'));
            BackgroundsParamII(z,c.ncube,c.fstar,c.vbc,c.vc,c.fx,c.sed,c.tau,c.zeta,c.feedback,c.delayParam,c.pop,c.fsfunc,~~c.photoHeating,c.photoHeating);
            rmpath('/a/home/cc/tree/taucc/students/physics/matanlotem/Work/SIM21/Clean_Fixed/');
            
            toc;
        end
        
        
        function deleteZ(z,c,dataPath,tmpDataPath)
            disp(['==Deleting z=',num2str(z),'==']);
            
            dataPath = p.dataPath;
            tmpDataPath = p.tmpDataPath;
            
            FdebugNaNs.deleteMagic(SIM21Analysis.xHIMagic,z,c,dataPath);
            FdebugNaNs.deleteMagic(SIM21Analysis.TKMagic,z-1,c,dataPath);
            FdebugNaNs.deleteMagic(SIM21Analysis.T21cmMagic,z,c,dataPath);
            FdebugNaNs.deleteMagic(SIM21Analysis.NeutMagic,z,c,dataPath);
            FdebugNaNs.deleteMagic(SIM21Analysis.EpsMagic,z,c,tmpDataPath);
            FdebugNaNs.deleteMagic(SIM21Analysis.XeMagic,z-1,c,tmpDataPath);
            FdebugNaNs.deleteMagic(SIM21Analysis.JLWMagic,z,c,tmpDataPath);
            FdebugNaNs.deleteMagic(SIM21Analysis.JalphaMagic,z,c,tmpDataPath);
            FdebugNaNs.deleteMagic(SIM21Analysis.LionMagic,z,c,tmpDataPath);
        end
        
        
        function deleteMagic(Magic,dataPath,c,z)
            fileName = SIM21Analysis.genDataFileName(dataPath,Magic,c.ID,z);
            if exist(fileName, 'file') == 2
                system(['rm ',fileName]);
            end
        end
        
        
        function compareZ(z,c,dataPath,tmpDataPath,compPath)
            disp(['==Checking z=',num2str(z),'==']);            
            FdebugNaNs.compareMagic(SIM21Analysis.xHIMagic,z,c,SIM21Analysis.xHIZ,dataPath,compPath);
            FdebugNaNs.compareMagic(SIM21Analysis.TKMagic,z-1,c,SIM21Analysis.TKZ,dataPath,compPath);
            FdebugNaNs.compareMagic(SIM21Analysis.T21cmMagic,z,c,SIM21Analysis.T21cmZ,dataPath,compPath);
            FdebugNaNs.compareMagic(SIM21Analysis.NeutMagic,z,c,SIM21Analysis.NeutZ,dataPath,compPath);
            FdebugNaNs.compareMagic(SIM21Analysis.EpsMagic,z,c,SIM21Analysis.EpsZ,tmpDataPath,compPath);
            FdebugNaNs.compareMagic(SIM21Analysis.XeMagic,z-1,c,SIM21Analysis.XeZ,tmpDataPath,compPath);
            FdebugNaNs.compareMagic(SIM21Analysis.JLWMagic,z,c,SIM21Analysis.JLWZ,tmpDataPath,compPath);
            FdebugNaNs.compareMagic(SIM21Analysis.JalphaMagic,z,c,SIM21Analysis.JalphaZ,tmpDataPath,compPath);
            FdebugNaNs.compareMagic(SIM21Analysis.LionMagic,z,c,SIM21Analysis.LionZ,tmpDataPath,compPath);
        end
        
        
        function compareMagic(Magic,z,c,zs,dataPath,compPath)
            if sum(find(zs==z))>0
                msg = '';
                a1 = SIM21Analysis.genDataFileName(dataPath,Magic,c.ID,z);
                a2 = SIM21Analysis.genDataFileName([dataPath,compPath],Magic,c.ID,z);
                if exist(a1, 'file') ~= 2
                    msg = [msg,'No File'];
                elseif exist(a2, 'file') ~= 2
                    msg = [msg,'No tmp File'];
                elseif isequal(importdata(a1),importdata(a2))
                    msg = [msg,'1'];
                else
                    msg = [msg,'0'];
                end
                disp(sprintf([Magic,'\t',msg,'\t',a1((length(dataPath)+1):end)]));
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Analayzing Stuff Tools
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function M = loadMat(Magic,zs,dataPath,c)
            disp(['==Loading ',Magic,'==']);
            M = zeros(max(zs),128,128,128);
            for z = zs
                fileName = SIM21Analysis.genDataFileName(dataPath,Magic,c.ID,z);
                if exist(fileName, 'file') == 2
                    M(z,:,:,:) = importdata(fileName);
                end
            end
        end
    end
end