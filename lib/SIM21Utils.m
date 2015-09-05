classdef SIM21Utils
    properties(Constant)
        libPath = '/a/home/cc/tree/taucc/students/physics/matanlotem/Work/SIM21/lib/';
        matrixPath = [SIM21Utils.libPath,'Matrices/'];

        %Data Matrix settings
        xHI = struct('z',[6:0.1:15,16:60],'magic','xHI','tmpData',0);
        TK = struct('z',[5:64],'magic','TK','tmpData',0);
        T21cm = struct('z',[6:60],'magic','T21cm','tmpData',0);
        Neut = struct('z',[6:0.1:15,16:60],'magic','Neut','tmpData',0);
        eps = struct('z',[6:60],'magic','eps','tmpData',1);
        xe = struct('z',[5:64],'magic','xe','tmpData',1);
        JLW = struct('z',[6,66],'magic','JLW','tmpData',1);
        Jalpha = struct('z',[6:60],'magic','Jalpha','tmpData',1);
        Lion = struct('z',[6:60],'magic','Lion','tmpData',1);

        
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


        function dataFileName = getDataFileName(c,dataType,z)
            % Get raw data matrix file name
            if ~ dataType.tmpData
                dataPath = c.dataPath;
            else
                dataPath = c.tmpDataPath;
            end
            dataFileName = [dataPath,dataType.magic,'_',num2str(z),c.ID,'.mat'];
        end


        function dataMat = getAllData(c,dataType)
            for z = 1:max(dataType.z)
                disp(num2str(z));
                dataMat(z,:,:,:) = zeros(128,128,128);
                if exist(SIM21Utils.getDataFileName(c,dataType,z)) == 2
                    dataMat(z,:,:,:) = importdata(SIM21Utils.getDataFileName(c,dataType,z));
                end
            end
        end


        function compareZ(c1,c2,z)
            disp(['==Checking z=',num2str(z),'==']);            
            SIM21Utils.compareMagic(c1,c2,z,SIM21Utils.xHI);
            SIM21Utils.compareMagic(c1,c2,z-1,SIM21Utils.TK);
            SIM21Utils.compareMagic(c1,c2,z,SIM21Utils.T21cm);
            SIM21Utils.compareMagic(c1,c2,z,SIM21Utils.Neut);
            SIM21Utils.compareMagic(c1,c2,z,SIM21Utils.eps);
            SIM21Utils.compareMagic(c1,c2,z-1,SIM21Utils.xe);
            SIM21Utils.compareMagic(c1,c2,z,SIM21Utils.JLW);
            SIM21Utils.compareMagic(c1,c2,z,SIM21Utils.Jalpha);
            SIM21Utils.compareMagic(c1,c2,z,SIM21Utils.Lion);
        end
        
        
        function compareMagic(c1,c2,z,dataType)
            if max(find(dataType.z==z))
                msg = '';
                a1 = SIM21Utils.getDataFileName(c1,dataType,z);
                a2 = SIM21Utils.getDataFileName(c2,dataType,z);
                if exist(a1, 'file') ~= 2
                    msg = [msg,'No File'];
                elseif exist(a2, 'file') ~= 2
                    msg = [msg,'No tmp File'];
                elseif isequal(importdata(a1),importdata(a2))
                    msg = [msg,'1'];
                else
                    msg = [msg,'0'];
                end
                disp(sprintf([dataType.magic,'\t',msg]));
            end
        end
    end
end