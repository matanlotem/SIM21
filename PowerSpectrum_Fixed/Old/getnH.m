% hydrogen proper number density
function nH = getnH(z)
    load Planck_parameters
    Y = 0.247; % Helium abundance by mass
    nH=(rhoc/mp)*(1-Y)*Ob*(1+z).^3; % 1/Mpc^3
end