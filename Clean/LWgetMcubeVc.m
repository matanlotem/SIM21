% Input parameters: 
% flag: with (1) or  without (0) vbc
% z: redshift
% JLW21 matrix of LW (delayed) intensities 
% flagM: here - circular velocity: 
% Vc = 4.2 or 3.7 km/s for molecular cooling, 
% Vc = 16.5 km/s atomic cooling
% Vc = 35 km/s heavy halos

function M = LWgetMcubeVc(flag, z, JLW21,Vcirc)

global vbc_cube
h = 0.6704;
Oc = 0.12038/h^2;
Ob = 0.022032/h^2;
Om = Ob+Oc;
Omz = Om*(1+z).^3./(Om.*(1+z).^3+1-Om);
d = Omz-1;
Dc = 18*pi^2+82*d-39*d^2;

a=4.015;
zrec = 1020;
Vc = (Vcirc^2 + (a*flag*vbc_cube*(1+z)/(1+zrec)*0.000097*3e5).^2).^0.5;

Mc = 10^8*(Vc/23.4).^3*((1+z)/10).^(-3/2)*(Om*Dc/Omz/18/pi^2).^-0.5/h;
M = Mc.*(1+ 6.96*(JLW21*4*pi).^0.47);



if (z<51)
    M = min(M,1e11*ones(size(M)));%Mass = logspace(log10(1e5),log10(7e8),NresM);
else
    M = min(M,1e8*ones(size(M)));
end

M = max(M,10^5*ones(size(M)));
Mc = [];   
 
end
