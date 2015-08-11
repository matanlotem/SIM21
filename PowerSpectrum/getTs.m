% calculate spin temperature
function Ts = getTs(TK, z0,ncube,fstar,flag,flagM,XeffTerm,Ispec,Reion,feedback,p,pop,FSfunc,photoheatingVersion)

Tcmb = 2.725*(1+z0); 
xA = getXA(TK, z0,ncube,fstar,flag,flagM,XeffTerm,Ispec,Reion,feedback,p,pop,FSfunc,photoheatingVersion);
xC= getXc(TK,z0); 
Tse = 0.4;%[K]

Ts = (1+xC+xA.*(1+Tse./TK).^(-1))./(Tcmb^(-1)+xC.*TK.^(-1)+xA.*TK.^(-1).*(1+Tse./TK).^(-1));


end