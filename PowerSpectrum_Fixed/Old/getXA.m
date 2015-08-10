function [x_talpha]=getXA(TK, z0,ncube,fstar,flag,flagM,XeffTerm,Ispec,Reion,feedback,p,pop,FSfunc,photoheatingVersion) 
%  z0=16;
%  ncube=0;
%  fstar=0.05; 
%  flag=0;
%  flagM=16.5;
%  XeffTerm=1;
%  Ispec=1;
%  Reion=0.075;
%  feedback=0;
%  p=0;
    load ('Planck_parameters');
    
    S_alpha2=getSalpha(TK, z0);
    JA_alpha = getLyA(TK, z0,ncube,fstar,flag,flagM,XeffTerm,Ispec,Reion,feedback,p,pop,FSfunc,photoheatingVersion);
    x_talpha=real(1.7e11*(1+z0).^(-1)*(S_alpha2.*JA_alpha));
end