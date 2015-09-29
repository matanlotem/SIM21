p = paramStudy();

for caseNum = [1,26:29]
    c = p.paramCases(caseNum).c;
    c.ncube = ncube;
    c.ID = c.getID();
end


xe = SIM21Utils.dataTypes.xe;
xHI = SIM21Utils.dataTypes.xHI;

if ~exist('xe29')
    xe29 = p.paramCases(29).c.getData(xe);
end
xeData = xe29;
if ~exist('xHI29')
    xHI29 = p.paramCases(29).c.getData(xHI);
end
xHIData = xHI29;

a1 = 64;
for z = xe.z(2:find(xe.z==max(xHI.z)))
    disp(z);
    outputName = [p.outputPath,'Graphs/xe/Case29_cube',num2str(ncube),'_a1_',num2str(a1),'_z_',num2str(z),'.png'];
    f = figure();   
    
    subplot(1,2,1,'Position',[0.025,0.2,0.45,0.55]);
    h = pcolor(1-squeeze(xeData(find(xe.z==z),a1,:,:)));
    h.EdgeColor = 'none';
    title(['xe - z = ',num2str(z)]);

    subplot(1,2,2,'Position',[0.525,0.2,0.45,0.55]);
    h = pcolor(squeeze(xHIData(find(xHI.z==z),a1,:,:)));
    h.EdgeColor = 'none';
    title(['xHI - z = ',num2str(z)]);

    saveas(f,outputName);
end