addpath('/a/home/cc/tree/taucc/students/physics/matanlotem/Work/SIM21/lib/');
p = paramStudy();

caseNum = 1;
c = p.paramCases(caseNum);
z = 62;

tmpDataPath = p.tmpDataPath;
dataPath = p.dataPath;
tmpDataPath = [p.tmpDataPath,'debug/'];
dataPath = [p.dataPath,'debug/'];

global pathname
global pathname_Data1
global pathname_Data2
global pathname_DataBackgrounds

pathname = '/scratch300/matanlotem/';
pathname_Data1 = tmpDataPath;
pathname_Data2 = dataPath;
pathname_DataBackgrounds = '/scratch300/matanlotem/DataBackgrounds_withPlanck/';

%FdebugNaNs.deleteZ(z,c,dataPath,tmpDataPath)
FdebugNaNs.runSimulation(c);
%FdebugNaNs.runZ(z,c)
for z=[64:-1:59]
    FdebugNaNs.compareZ(z,c,dataPath,tmpDataPath);
end

%disp('==Loading TK==')
%TKData = zeros(64,128,128,128);
%for z = 5:64
%    TKData(z,:,:,:) = importdata(SIM21Analysis.genDataFileName(p.dataPath,'TK',p.paramCases(caseNum).ID,z));
%end
%disp('==Done Loading==');
%TKRows = TKData(:,:);
%
%cat(1,mean(TKRows,2)',[(sum(diff(TKRows)<0,2)/128^3)',[0]])
