classdef SIM21Utils
    methods(Static)
        function M = importMatrix(MName)
            M = importdata(SIM21Utils.matrixPath(MName));
        end
        
        
        function Mpath = matrixPath(MName)
            Mpath=['Matrices/',MName,'.mat'];
        end
        
        
        function ID = getID(MyCube,MyStar,MyVBC,MyVc,MyFX,MySED,MyTau,MyFeed,DelayParam,MyPop,FSfunc,photoheatingVersion)
            ID = ['_' num2str(MyCube)...
                  '_' num2str(MyStar) '_' num2str(MyVBC) '_' num2str(MyVc)...
                  '_' num2str(MyFX) '_' num2str(MySED) '_' num2str(MyTau)...
                  '_' num2str(MyFeed) '_' num2str(DelayParam) '_' num2str(MyPop) '_' num2str(FSfunc) '_' num2str(photoheatingVersion)];
        end
    end
end