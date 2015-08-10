function tGP = getTauGPn(z) 
    load Planck_parameters
    gamma = 50*10^6; % 1/s
    nH = getnH(z);% 1/Mpc^3
    tGP=(3*nH*(lambda_a)^3*gamma)./(2*getHz(z)*3.241e-20);
end
