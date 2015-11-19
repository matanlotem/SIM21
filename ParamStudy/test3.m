res = [];
for caseNum = 50:65
	disp(caseNum)
	runCase = p.paramCases(caseNum).c;
	caseTK = runCase.getData('TK');
	caseNeut = runCase.getData('Neut');
	[z,r,c,d] = ind2sub(size(caseTK(1:40,:,:,:)),find(isnan(caseTK(1:40,:,:,:))));
	nind = SIM21Utils.matchDataZIndexes(z,'TK','Neut');
	a = [nind',r,c,d];
	b = a(find(nind~=0),:,:,:);
	g = isempty(find(caseNeut(sub2ind(size(caseNeut),b(:,1),b(:,2),b(:,3),b(:,4)))~=0));
	res = [res;[caseNum,g,length(z)]];
	disp(res(end,:));
end