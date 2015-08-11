% Hubble constant
function Hz = getHz(z)
% in km/s/Mpc
load Planck_parameters

Hz = H0*sqrt(Om*(1+z).^3+OLambda+8.5522e-05*(1+z).^4); 

end