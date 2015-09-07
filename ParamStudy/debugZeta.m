p = paramStudy();

for caseNum = [34,35,42,43]
	c = p.paramCases(caseNum).c;
	c.setPathExt('Test8');
	SIM21Analysis.plotGraphsByZ([p.paramCases(1).c,c]);
	SIM21Analysis.calcSpecialParams(c);
end