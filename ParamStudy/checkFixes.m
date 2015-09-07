p = paramStudy();
cases(1) = p.paramCases(1).c;
cases(2) = copy(cases(1));
cases(2).setPathExt('Test2');
cases(2).name = 'Case 1 - 3 Fixes';
cases(2).setOutputPath('/scratch300/matanlotem/ParamStudy/Tests/Case_1_3_Fixes/',1);
cases(3) = copy(cases(1));
cases(3).setPathExt('Test4');
cases(3).name = 'Case 1 - All Fixes';
cases(3).setOutputPath('/scratch300/matanlotem/ParamStudy/Tests/Case_1_All_Fixes/',1);
%SIM21Analysis.plotGraphsByZ(cases);

cases50(1) = p.paramCases(1).c;
cases50(2) = p.paramCases(50).c;
cases50(3) = copy(cases50(2));
cases50(3).name = 'Case 50 - All Fixes';
cases50(3).setPathExt('Test5');
cases50(3).setOutputPath('/scratch300/matanlotem/ParamStudy/Tests/Case_50_All_Fixes/',1);
%SIM21Analysis.plotGraphsByZ(cases50);
%for i=[2:3]
%	SIM21Analysis.calcSpecialParams(cases50(i));
%end

cases51(1) = p.paramCases(1).c;
cases51(2) = p.paramCases(51).c;
cases51(3) = copy(cases51(2));
cases51(3).name = 'Case 51 - All Fixes';
cases51(3).setPathExt('Test5');
cases51(3).setOutputPath('/scratch300/matanlotem/ParamStudy/Tests/Case_51_All_Fixes/',1);
%SIM21Analysis.plotGraphsByZ(cases51);
%for i=[2:3]
%	SIM21Analysis.calcSpecialParams(cases51(i));
%end

cases10(1) = p.paramCases(1).c;
cases10(2) = p.paramCases(10).c;
cases10(3) = copy(cases10(2));
cases10(3).name = 'Case 10 - 3 Fixes';
cases10(3).setPathExt('Test3');
cases10(3).setOutputPath('/scratch300/matanlotem/ParamStudy/Tests/Case_10_3_Fixes/',1);
cases10(4) = copy(cases10(2));
cases10(4).name = 'Case 10 - All Fixes';
cases10(4).setPathExt('Test6');
cases10(4).setOutputPath('/scratch300/matanlotem/ParamStudy/Tests/Case_10_All_Fixes/',1);
%for i=[2:4]
%	SIM21Analysis.calcSpecialParams(cases10(i));
%end
%SIM21Analysis.plotGraphsByZ(cases10);

cases20(1) = p.paramCases(1).c;
cases20(2) = p.paramCases(20).c;
cases20(3) = copy(cases20(2));
cases20(3).name = 'Case 20 - 3 Fixes';
cases20(3).setPathExt('Test3');
cases20(3).setOutputPath('/scratch300/matanlotem/ParamStudy/Tests/Case_20_3_Fixes/',1);
cases20(4) = copy(cases20(2));
cases20(4).name = 'Case 20 - All Fixes';
cases20(4).setPathExt('Test6');
cases20(4).setOutputPath('/scratch300/matanlotem/ParamStudy/Tests/Case_20_All_Fixes/',1);
%for i=[2:4]
%	SIM21Analysis.calcSpecialParams(cases20(i));
%end
%SIM21Analysis.plotGraphsByZ(cases20);

cases30(1) = p.paramCases(1).c;
cases30(2) = p.paramCases(30).c;
cases30(3) = copy(cases30(2));
cases30(3).name = 'Case 30 - 3 Fixes';
cases30(3).setPathExt('Test3');
cases30(3).setOutputPath('/scratch300/matanlotem/ParamStudy/Tests/Case_30_3_Fixes/',1);
cases30(4) = copy(cases30(2));
cases30(4).name = 'Case 30 - All Fixes';
cases30(4).setPathExt('Test6');
cases30(4).setOutputPath('/scratch300/matanlotem/ParamStudy/Tests/Case_30_All_Fixes/',1);
%for i=[2:4]
%	SIM21Analysis.calcSpecialParams(cases30(i));
%end
%SIM21Analysis.plotGraphsByZ(cases30);

cases60(1) = p.paramCases(1).c;
cases60(2) = p.paramCases(60).c;
cases60(3) = copy(cases60(2));
cases60(3).name = 'Case 60 - 3 Fixes';
cases60(3).setPathExt('Test3');
cases60(3).setOutputPath('/scratch300/matanlotem/ParamStudy/Tests/Case_60_3_Fixes/',1);
cases60(4) = copy(cases60(2));
cases60(4).name = 'Case 60 - All Fixes';
cases60(4).setPathExt('Test6');
cases60(4).setOutputPath('/scratch300/matanlotem/ParamStudy/Tests/Case_60_All_Fixes/',1);
%for i=[2:4]
%	SIM21Analysis.calcSpecialParams(cases60(i));
%end
%SIM21Analysis.plotGraphsByZ(cases60);

cases55(1) = p.paramCases(1).c;
cases55(2) = p.paramCases(55).c;
cases55(3) = copy(cases55(2));
cases55(3).name = 'Case 55 - All Fixes';
cases55(3).setPathExt('Test7');
cases55(3).setOutputPath('/scratch300/matanlotem/ParamStudy/Tests/Case_55_All_Fixes/',1);
%for i=[2:3]
%	SIM21Analysis.calcSpecialParams(cases55(i));
%end
%SIM21Analysis.plotGraphsByZ(cases55);

cases59(1) = p.paramCases(1).c;
cases59(2) = p.paramCases(59).c;
cases59(3) = copy(cases59(2));
cases59(3).name = 'Case 59 - All Fixes';
cases59(3).setPathExt('Test7');
cases59(3).setOutputPath('/scratch300/matanlotem/ParamStudy/Tests/Case_59_All_Fixes/',1);
%for i=[2:3]
%	SIM21Analysis.calcSpecialParams(cases59(i));
%end
%SIM21Analysis.plotGraphsByZ(cases59);

%%% check xe-XHI relations
for c = [cases50(3),cases51(3),cases55(3),cases59(3),cases60(4)];
	figSettings = SIM21Analysis.initFigSettings('xe & xHI (z)',c.name,'z + 1','x');
	figSettings.lines{end+1} = SIM21Analysis.plotLine(SIM21Analysis.zPlus1(SIM21Analysis.getZData(c,'xHI')),'xHI');
	xe1 = SIM21Analysis.getZData(c,'xe');
	xe1(2,:) = 1-xe1(2,:);
	figSettings.lines{end+1} = SIM21Analysis.plotLine(SIM21Analysis.zPlus1(xe1),'1-xe');
	figName = [c.outputPath,'xe-xHI',c.ID,'.png'];
	SIM21Analysis.plotData(figName,figSettings);
end