classdef SIM21Testing
    properties(Constant)
	end

	methods(Static)
		function Test(a,b,c)
			if exist('a','var')
				disp(a)
			end
			if exist('b','var')
				disp(b)
			end
			if exist('c','var')
				disp(c)
			end
		end

		function Test1(varargin)
			SIM21Testing.Test(varargin{:});
		end

		function Test2(a,b)
			disp(varargin)
		end


		function compareZ(c1,c2,z)
            disp(['==Checking z=',num2str(z),'==']);            
            SIM21Testing.compareMagic(c1,c2,z,SIM21Utils.dataTypes.xHI);
            SIM21Testing.compareMagic(c1,c2,z-1,SIM21Utils.dataTypes.TK);
            SIM21Testing.compareMagic(c1,c2,z,SIM21Utils.dataTypes.T21cm);
            SIM21Testing.compareMagic(c1,c2,z,SIM21Utils.dataTypes.Neut);
            SIM21Testing.compareMagic(c1,c2,z,SIM21Utils.dataTypes.eps);
            SIM21Testing.compareMagic(c1,c2,z-1,SIM21Utils.dataTypes.xe);
            SIM21Testing.compareMagic(c1,c2,z,SIM21Utils.dataTypes.JLW);
            SIM21Testing.compareMagic(c1,c2,z,SIM21Utils.dataTypes.Jalpha);
            SIM21Testing.compareMagic(c1,c2,z,SIM21Utils.dataTypes.Lion);
        end
        
        
        function compareMagic(c1,c2,z,dataType)
        	dataType = SIM21Utils.getDataType(dataType);
            if max(find(dataType.z==z))
                msg = '';
                fileName1 = SIM21Utils.getDataFileName(c1,dataType,z);
                fileName2 = SIM21Utils.getDataFileName(c2,dataType,z);
                if exist(fileName1, 'file') ~= 2
                    msg = [msg,'No File 1 ',fileName1];
                elseif exist(fileName2, 'file') ~= 2
                    msg = [msg,'No File 2',fileName2];
                else
                	data1 = importdata(fileName1);
                	data2 = importdata(fileName2);
                	
                	diffs1 = data1(data1 ~= data2);
                	diffs2 = data2(data1 ~= data2);
                	if isequal(data1,data2)
                    	msg = [msg,'1'];
                    elseif ~max(abs((diffs1-diffs2) ./ (diffs1+diffs2)) > 10^-6)
                   		msg = [msg,'2'];
                	else
                    	msg = [msg,'0'];
                    end
                end
                disp(sprintf([dataType.magic,'\t',msg]));
            end
        end
	end
end