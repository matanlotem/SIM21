function Rc = getHorizon(zcenter)
 
load Planck_parameters
%  N=100000000;
% z=linspace(zcenter,10^8,N);
% H=H0*sqrt(Om*(1+z).^3+OLambda);
% dz = c./H; 
% Rc = trapz(z,dz); %Mpc

F = @(x)c./(H0*sqrt(Om*(1+x).^3+OLambda));
%Rc = integral(F,zcenter,Inf);
Rc = quadl(F,zcenter,1e10);
end