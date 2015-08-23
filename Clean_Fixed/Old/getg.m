function g = getg(z,A,B,Lion,eps)
   fHe=6/(6+76);
   kb = (1.38e-16/1.6e-12); % eV/K
   dtdz =-3.09e19./(getHz(z).*(1+z)); % in seconds
   g = 2-getf(z,A,B,Lion)+(2/3/kb)*(1+z).*dtdz*eps./(exp(B).*exp(A))+((exp(A)-1)./(exp(A)+fHe)).*8.55e-13.*(2.725*(1+z)-exp(B)).*(1+z).^5./(exp(B));

end