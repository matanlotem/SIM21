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


		function flag = compareZ(c1,c2,z)
            disp(['==Checking z=',num2str(z),'==']);
            flag = 1;
            flag = flag && SIM21Testing.compareMagic(c1,c2,z,SIM21Utils.dataTypes.xHI);
            flag = flag && SIM21Testing.compareMagic(c1,c2,z-1,SIM21Utils.dataTypes.TK);
            flag = flag && SIM21Testing.compareMagic(c1,c2,z,SIM21Utils.dataTypes.T21cm);
            flag = flag && SIM21Testing.compareMagic(c1,c2,z,SIM21Utils.dataTypes.Neut);
            flag = flag && SIM21Testing.compareMagic(c1,c2,z,SIM21Utils.dataTypes.eps);
            flag = flag && SIM21Testing.compareMagic(c1,c2,z-1,SIM21Utils.dataTypes.xe);
            flag = flag && SIM21Testing.compareMagic(c1,c2,z,SIM21Utils.dataTypes.JLW);
            flag = flag && SIM21Testing.compareMagic(c1,c2,z,SIM21Utils.dataTypes.Jalpha);
            flag = flag && SIM21Testing.compareMagic(c1,c2,z,SIM21Utils.dataTypes.Lion);
        end
        
        
        function flag = compareMagic(c1,c2,z,dataType)
        	dataType = SIM21Utils.getDataType(dataType);
            flag = 0;
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
                        flag = 1;
                    elseif ~max(abs((diffs1-diffs2) ./ (diffs1+diffs2)) > 10^-6)
                        flag = 2;
                    end
                    msg = [msg,num2str(flag)];
                end
                disp(sprintf([dataType.magic,'\t',msg]));
            else
                flag=3;
            end
        end


        function testSim(c,caseName,zs)
            runCase = copy(c);
            sampleCase = copy(c);
            runCase.setPathExt(['TESTING/WORK/',caseName,'/']);
            sampleCase.setPathExt(['TESTING/SAMPLES/',caseName,'/']);

            if ~exist('zs','var')
                zs = [6,12.5,15,30,60];
            end

            res = zeros(1,length(zs));
            timings = zeros(1,length(zs));
            for z = zs
                tic;
                SIM21Utils.runZ(runCase,z);
                res(zs==z) = SIM21Testing.compareZ(runCase,sampleCase,z);
                timings(zs==z) = toc;
            end
            disp(cat(1,zs,res,timings/60));
        end
	end
end