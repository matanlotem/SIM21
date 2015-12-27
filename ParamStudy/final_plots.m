p = paramStudy();
regularCase = p.regularCase;
bigZetaCases = [58:65];
metalCoolingCases = [37:2:49];
badTSCases = setdiff([],[bigZetaCases,metalCoolingCases]);
smallVarCases = setdiff(p.smallVarCases,[bigZetaCases,metalCoolingCases]);
largeVarCases = setdiff(p.largeVarCases,[bigZetaCases,metalCoolingCases]);

figSettings = SIM21Analysis.initFigSettings('','','T21cm [mK]','1+\ni');
