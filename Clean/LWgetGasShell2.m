function fgas_shell = LWgetGasShell2(fgas_wk,Rmin,Rmax)

%%% smooth fgas field 
if(Rmin > 0)
    f1 = LWsmooth2(fgas_wk,Rmin); %gives average mass per pixel in stars
    f2 = LWsmooth2(fgas_wk,Rmax);
    Vol1 = 4/3*pi*Rmin^3;
    Vol2 = 4/3*pi*Rmax^3;
    fgas_shell = (f2*Vol2 - f1*Vol1); % gives total mass in stars in shell in solar mass
else
    f2 = LWsmooth2(fgas_wk,Rmax);
    Vol2 = 4/3*pi*Rmax^3;
    fgas_shell = f2*Vol2; % gives total mass in stars in shell in solar mass
end