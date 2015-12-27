%p = paramStudy();
%outputPath = [p.outputPath,'Graphs/'];
%for perc = [25,50,75]
%	xfield = {'specialParams',['xHI',num2str(perc)],'z'};
%	%yfield = {'tempParams','minMaxSlope','zDiff'};
%	yfield = {'paramCases','atau'};
%	xlabel = strjoin(xfield(2:end),'.');
%	ylabel = strjoin(yfield(2:end),'.');
%	figSettings = SIM21Analysis.initFigSettings([ylabel,'(',xlabel,')'],'',xlabel,ylabel);
%	caseNums = p.workCases;
%	figSettings.plots = [figSettings.plots, p.getLines(xfield,yfield,{'paramCases','c','vc'},caseNums)];
%	%figSettings.plots = [figSettings.plots, p.getLines(xfield,yfield,{'paramCases','c','fstar'},caseNums)];
%	figSettings.xTick = [1:30];
%	p.plotSigScatter(xfield,yfield,figSettings,[outputPath,ylabel,'_',xlabel,'.png'],caseNums);
%end

%%yfield = {'specialParams','maxSlope','z'};
%%xfield = {'paramCases','c','fstar'};
%%xlabel = strjoin(xfield(2:end),'.');
%%ylabel = strjoin(yfield(2:end),'.');
%%figSettings = SIM21Analysis.initFigSettings([ylabel,'(',xlabel,')'],'',xlabel,ylabel);
%%figSettings.log = 'x';
%%figSettings.xTick = [0.005,0.015,0.05,0.158,0.5];
%%caseNums = [1,66:69];
%%p.plotSigScatter(xfield,yfield,figSettings,[outputPath,ylabel,'_',xlabel,'_66-69.png'],caseNums);
%%
%%yfield = {'specialParams','maxSlope','z'};
%%xfield = {'paramCases','c','vc'};
%%xlabel = strjoin(xfield(2:end),'.');
%%ylabel = strjoin(yfield(2:end),'.');
%%figSettings = SIM21Analysis.initFigSettings([ylabel,'(',xlabel,')'],'',xlabel,ylabel);
%%figSettings.log = 'y';
%%%figSettings.xTick = [0.005,0.015,0.05,0.158,0.5];
%%caseNums = [70:73];
%%p.plotSigScatter(xfield,yfield,figSettings,[outputPath,ylabel,'_',xlabel,'_70-73.png'],caseNums);

p = paramStudy();
outputPath = [p.outputPath,'Graphs/zT/'];
caseNums = [p.regularCase,p.smallVarCases,p.largeVarCases];
xfield = {'specialParams','minT21cm','z'};
yfield = {'specialParams','minT21cm','T'};
prefix = 'c';
[figSettings,outputName] = p.initSSFigSettings(xfield,yfield,outputPath,[prefix,'_vc_']);
figSettings.plots = [figSettings.plots, p.getLines(xfield,yfield,{'paramCases','c','vc'},caseNums)];
p.plotSigScatter(xfield,yfield,figSettings,outputName,caseNums)
[figSettings,outputName] = p.initSSFigSettings(xfield,yfield,outputPath,[prefix,'_fstar_']);
figSettings.plots = [figSettings.plots, p.getLines(xfield,yfield,{'paramCases','c','fstar'},caseNums)];
p.plotSigScatter(xfield,yfield,figSettings,outputName,caseNums)
[figSettings,outputName] = p.initSSFigSettings(xfield,yfield,outputPath,[prefix,'_fx_']);
figSettings.plots = [figSettings.plots, p.getLines(xfield,yfield,{'paramCases','c','fx'},caseNums)];
p.plotSigScatter(xfield,yfield,figSettings,outputName,caseNums)
[figSettings,outputName] = p.initSSFigSettings(xfield,yfield,outputPath,[prefix,'_sed_']);
figSettings.plots = [figSettings.plots, p.getLines(xfield,yfield,{'paramCases','c','sed'},caseNums)];
p.plotSigScatter(xfield,yfield,figSettings,outputName,caseNums)
[figSettings,outputName] = p.initSSFigSettings(xfield,yfield,outputPath,[prefix,'_tau_']);
figSettings.plots = [figSettings.plots, p.getLines(xfield,yfield,{'paramCases','c','tau'},caseNums)];
p.plotSigScatter(xfield,yfield,figSettings,outputName,caseNums)
%p.plotBasicSigScatter({'specialParams','minSlope','z'},{'specialParams','minSlope','T'},outputPath);
%p.plotBasicSigScatter({'specialParams','maxSlope','z'},{'specialParams','maxSlope','T'},outputPath);
%p.plotBasicSigScatter({'specialParams','minT21cm','z'},{'specialParams','minT21cm','T'},outputPath);
%p.plotBasicSigScatter({'specialParams','maxT21cm','z'},{'specialParams','maxT21cm','T'},outputPath);