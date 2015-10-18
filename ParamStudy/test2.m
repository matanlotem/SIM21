addpath('/a/home/cc/tree/taucc/students/physics/matanlotem/Work/SIM21/lib/');
p = paramStudy();
pt = paramStudy();
for caseNum = [1:73]
	c = pt.paramCases(caseNum).c;
	runName = ['PS_New_',num2str(caseNum)];
	c.setPathExt([runName,'_New']);
	c.setOutputPath([c.outputPath(1:end-1),'_New/']);
end
ptorkCases = pt.areRun([pt.regularCase,pt.smallVarCases,pt.largeVarCases,pt.otherCases]);

ct1 = pt.paramCases(1).c;
for caseNum = [1:9]
	c = p.paramCases(caseNum).c;
	ct = pt.paramCases(caseNum).c;
	SIM21Analysis.plotAllGraphsByZ([c,ct1,ct]);
end

%pt.plot4CasesZGraphs();