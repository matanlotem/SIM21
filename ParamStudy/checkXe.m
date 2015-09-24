p = paramStudy();
xe = SIM21Utils.dataTypes.xe;
xHI = SIM21Utils.dataTypes.xHI;
if ~exist('xe35')
	xe35 = p.paramCases(35).c.getData(xe);
end
if ~exist('xHI35')
	xHI35 = p.paramCases(35).c.getData(xHI);
end
a1 = 64;
for z = xe.z(2:end)
	disp(z);
	outputName = [p.outputPath,'Graphs/xe/a1_',num2str(a1),'_z_',num2str(z),'.png'];
	f = figure();	
	
	subplot(1,2,1,'Position',[0.025,0.2,0.45,0.55]);
	h = pcolor(1-squeeze(xe35(find(xe.z==z),a1,:,:)));
	h.EdgeColor = 'none';
	title(['xe - z = ',num2str(z)]);

	subplot(1,2,2,'Position',[0.525,0.2,0.45,0.55]);
	h = pcolor(squeeze(xHI35(find(xHI.z==z),a1,:,:)));
	h.EdgeColor = 'none';
	title(['xHI - z = ',num2str(z)]);

	saveas(f,outputName);
end