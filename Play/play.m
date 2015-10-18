tic;
n = 64;
R = 64;
box = rand(n);
box = ones(n);

combox = repmat(box(),[1,1,n,n,R]);

[x1,y1,x2,y2,r] = ndgrid(1:n,1:n,1:n,1:n,1:R);
rind = ceil(((x1-x2).^2 + (y1-y2).^2).^0.5) == r;
a = squeeze(sum(sum(sum(combox.*rind,3),4),5));
toc;