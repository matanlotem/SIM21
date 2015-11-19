p = paramStudy();
f = findZeta();
for caseNum = [1]
	caseNum
	c = p.paramCases(caseNum).c;
	runName = ['PS_',num2str(caseNum)];
	c.setPathExt([runName]);
	c.ID = c.newID();
	SIM21Utils.runSimulation(c,runName);
	%zetaCase = f.getZetaCase(caseNum,c.zeta);
	%f.copyZetaCase(zetaCase,c);
	%SIM21Utils.changeID(c,c.newID());
end
%caseNum = 39;
%c = p.paramCases(caseNum).c;
%runName = ['PS_',num2str(caseNum)];
%c.setPathExt(runName);
%SIM21Utils.copyCaseData(f.getZetaCase(38,26).c,c,1);
%c1 = copy(c);
%c1.tau = 0.066;
%c1.ID = c1.getID();
%SIM21Utils.changeID(c1,c.ID);