classdef SIM21Utils
    properties(Constant)
        libPath = '/a/home/cc/tree/taucc/students/physics/matanlotem/Work/SIM21/lib/';
        matrixPath = [SIM21Utils.libPath,'Matrices/'];

        %Data Matrix settings
        xHI = struct('z',[6:0.1:15,16:60],'magic','xHI');
        %xHIMagic = 'xHI';
        %xHIZ = [6:0.1:15,16:60];
        TK = struct('z',[5:64],'magic','TK');
        %TKMagic = 'TK';
        %TKZ = [5:64];
        T21cm = struct('z',[6:60],'magic','T21cm');
        %T21cmMagic = 'T21cm';
        %T21cmZ = [6:60];
        Neut = struct('z',[6:0.1:15,16:60],'magic','Neut');
        %NeutMagic = 'Neut';
        %NeutZ = [6:0.1:15,16:60];
        eps = struct('z',[6:60],'magic','eps');
        %EpsMagic = 'eps';
        %EpsZ = [6:60];
        xe = struct('z',[5:64],'magic','xe');
        %XeMagic = 'xe';
        %XeZ = [5:64];
        JLW = struct('z',[6,66],'magic','JLW');
        %JLWMagic = 'JLW';
        %JLWZ = [6:66];
        Jalpha = struct('z',[6:60],'magic','Jalpha');
        %JalphaMagic = 'Jalpha';
        %JalphaZ = [6:60];
        Lion = struct('z',[6:60],'magic','Lion');
        %LionMagic = 'Lion';
        %LionZ = [6:60];

        
    end
    
    methods(Static)
        function M = importMatrix(MName)
            M = importdata(SIM21Utils.getMatrixPath(MName));
        end
        
        
        function Mpath = getMatrixPath(MName)
            Mpath=[SIM21Utils.matrixPath,MName,'.mat'];
        end
        
        
        function ID = getID(MyCube,MyStar,MyVBC,MyVc,MyFX,MySED,MyTau,MyFeed,DelayParam,MyPop,FSfunc,photoheatingVersion)
            ID = ['_' num2str(MyCube)...
                  '_' num2str(MyStar) '_' num2str(MyVBC) '_' num2str(MyVc)...
                  '_' num2str(MyFX) '_' num2str(MySED) '_' num2str(MyTau)...
                  '_' num2str(MyFeed) '_' num2str(DelayParam) '_' num2str(MyPop) '_' num2str(FSfunc) '_' num2str(photoheatingVersion)];
        end
    end
end