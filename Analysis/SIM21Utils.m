classdef SIM21Utils
    methods(Static)
        function M = importMatrix(MName)
            M = importdata(['Matrices/',MName,'.mat']);
        end
    end
end