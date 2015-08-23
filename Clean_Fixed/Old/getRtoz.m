% output: z+dz
function zs = getRtoz(R,z)

H0=67.04;
c=3e5;
omm=.3169;
k = 2*c/(H0*sqrt(omm));
Rcrit = k./sqrt(1+z);

%if (R<min(Rcrit,getHorizon(z)))
if (R<getHorizon(z))
    zs = (k^2*z + 2*k*R*sqrt(1+z) - R.^2.*(1+z))./(k^2 - 2*k*R*sqrt(1+z) + R.^2*(1+z));
else 
    zs = 1e15;
end
% zs = (k^2*z + 2*k*R*sqrt(1+z) - R.^2.*(1+z))./(k^2 - 2*k*R*sqrt(1+z) + R.^2*(1+z)).*(R<Rcrit)+...
%    (k^2*z + 2*k*R(IndR)*sqrt(1+z) - R(IndR).^2.*(1+z))./(k^2 - 2*k*R(IndR)*sqrt(1+z) + R(IndR).^2*(1+z)).*(R>Rcrit);

end