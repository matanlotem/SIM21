function f = getf(z,A,B,Lion)
   
    global delta_cube
    
    
    load Planck_parameters
    
    dtdz = -3.09e19./(getHz(z).*(1+z)); % in seconds
    fHe = 1-0.921;
    fHI = 0.921;
    nb0 = 2.1e-7*(1+fHe/fHI);
    alphaA=4.2e-13*(exp(B)/10^4).^(-0.76);
    f = (1+z).*dtdz*Lion./exp(A)-(1+z)^4.*dtdz.*fHI*nb0.*alphaA.*((exp(A)-1).^2./exp(A)).*(1+delta_cube*LWgetDz(z)/D40);

end