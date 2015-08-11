classdef SIM21Utils
    methods(Static)
        function M = importMatrix(MName)
            M = importdata(SIM21Utils.matrixPath(MName));
        end
        
        
        function Mpath = matrixPath(MName)
            Mpath=['Matrices/',MName,'.mat'];
        end
    end
end