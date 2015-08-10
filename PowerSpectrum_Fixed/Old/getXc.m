% collisional coupling coefficient
function Xc = getXc(TK,z0)
    global delta_cube
    load Planck_parameters
    nH = (ones(size(delta_cube))+max(-1,min(1,delta_cube))).*getnH(z0)/(3.1e24)^3;% 1/cm^3
    kappa = getKappa(TK);%cm^3/s
    Tcmb = 2.725*(1+z0);%K
    Tstar = 0.068; %K
    Xc = (nH.*kappa/A10)*(Tstar/Tcmb);
end