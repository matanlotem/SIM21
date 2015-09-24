p = paramStudy();
for caseNum = [1,26:29]
	c = p.paramCases(caseNum).c;
	c.ncube = 2;
	c.ID = c.getID();
	runName = ['paramStudy_',num2str(caseNum)];
	c.setPathExt(runName);
	SIM21Utils.runSimulation(c,runName);
end