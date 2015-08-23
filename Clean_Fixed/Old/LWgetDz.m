function y = LWgetDz(z)  % calculate growth function

a=1./(1+z);

global omm 
omm=.3169; % omega matter
oml=1-omm;
%ap=linspace(0,a,100);
%y = sqrt(oml*a^3 + omm)/a^(3/2)*trapz(ap,ap.^(3/2)./(oml*ap.^3+omm).^(3/2));
y = sqrt(oml*a^3 + omm)/a^(3/2)*quad(@IntegrandDz,0,a,1e-4);
%normalize such that D(0) = 1.
y=y*1.004;%/1.125940274656245;

end

