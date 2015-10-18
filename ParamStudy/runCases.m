p = paramStudy();
for caseNum = [6:9]
	c = p.paramCases(caseNum).c;
	%c.ncube = 1;
	%c.ID = c.getID();
	runName = ['PS_New_',num2str(caseNum)];
	c.setPathExt([runName,'_New']);
	SIM21Utils.runSimulation(c,runName);
end