%
%% Calculates the critical mass for the photoheating effect (version 2)
%% + directory - Where the xHI files are stored.
%% + filePrefix - A.K.A. '_xHI_'.
%% + zMax - Maximal value of z in the simulation for which the function
%%          will be called upon.
%% + zMin - Minimal value of z in the simulation for which the function
%%          will be called upon.
%% + N - Mass matrix dimension size.
%% + z - The redshift to calculate for.
%%
%% Returns:
%% + Mcrit - NxNxN matrix of the critical glaxy mass
%% + xHI - 1-ionization (as a fraction) for z+1
%% + currentIonization - Ionization (as a fraction) for z+1
%%
%% To run the tests set DEBUG=1 (define globally)
function [Mcrit,currentIonization,a,b] = calculateMcritV2(flag,flagM,JLW21,fstarM,fstarA,FSfunc,directory, filePrefix, zMax, zMin, N, z)

global DEBUG;
if isempty(DEBUG)
    DEBUG = 0;
end
global gridInterpKey;

% persistent zer;
% if isempty(zer)
    zer = zeros(N,N,N);
% end

Rset = [linspace(1e-10,500/3,110) logspace(log10(550/3),log10(15000/3),20)];
Rmin = Rset(1:end-1);
Rmax = Rset(2:end);
Nshells = length(Rmin);
Lpix=3;

ZN = zMax-15+((15-zMin)*10+1); % Number of z values - this should be globally defined!!

zs=[60:-1:16 15:-0.1:6];
previousZIndex=find((round(zs*10)/10)==(round(z*10)/10))-1;



    load([directory filePrefix num2str(zs(previousZIndex-1)) gridInterpKey]);
    z2xHI=xHI;
    load([directory filePrefix num2str(zs(previousZIndex)) gridInterpKey]);
    z1xHI=xHI;
    x_interp(1,:,:,:)=z1xHI;
    x_interp(2,:,:,:)=z2xHI;
    currentxInter=(squeeze(10.^interp1(log10([zs(previousZIndex),zs(previousZIndex-1)]),log10(x_interp+1e-10),log10(z),'linear','extrap')));
    currentxInter=currentxInter.*logical(squeeze(x_interp(1,:,:,:))>0);
    x_interp=[];
    
   

    xHI1=max(0,min(1,(currentxInter).*0.9)); 
    xHI2=max(0,min(1,(currentxInter).*1.1)); 
    

    halfCurrentIonization=(1-xHI2)/2;
    MinHalf=min(min(min(halfCurrentIonization)));



i=0;
mark=0;
while mark==0 && i~=previousZIndex
    load([directory filePrefix num2str(zs(previousZIndex-i)) gridInterpKey]);
    if max(max(max(1-xHI)))<MinHalf
        mark=i;
    end
    i=i+1;
end

if mark==0
    mark=previousZIndex-1;
end

ionizationHistory=zeros(N,N,N,mark+1);
ind=1;
for i=(previousZIndex-mark):previousZIndex
        load([directory filePrefix num2str(zs(i)) gridInterpKey]);
        ionizationHistory(:,:,:,ind) = 1-xHI;
        ind=ind+1;
end

fin=mark+1;







if z<(zMax-1)

    %%
    %% Initialize
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%
    %% Find zIN and zIN_ion
    %% + To optimize only compare between the previous 1/2 ionization
    %%   index and the current z-1 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
   
    
    %% Load current ionization



  
    Mcrit=calcMcritV2Fast(z,zMax,N,zs,1-xHI1,fin, ionizationHistory)  ;
    
    
    fgas1=fgasForPhotoV2(flag,flagM,JLW21,z,fstarM,fstarA,FSfunc,1-xHI1,Mcrit)/fstarM;
  
    currentxInter=[];
        
    Mcrit=calcMcritV2Fast(z,zMax,N,zs,1-xHI2,fin, ionizationHistory)  ;

    
    fgas2=fgasForPhotoV2(flag,flagM,JLW21,z,fstarM,fstarA,FSfunc,1-xHI2,Mcrit)/fstarM;
    
    eq=(xHI1==xHI2);
    
    a = zeros(N,N,N);
    b = zeros(N,N,N);
    a(eq) = 0;
    b(eq) = fgas1(eq);
    a(~eq) = (fgas2(~eq)-fgas1(~eq))./(xHI2(~eq)-xHI1(~eq));
    b(~eq) = (xHI1(~eq).*fgas2(~eq)-xHI2(~eq).*fgas1(~eq))./(xHI1(~eq)-xHI2(~eq));
    
    zeta=fstarM/0.05*19.48;
    
    xHI = min(1,max(0,(1-b.*zeta)./(1+a.*zeta))); % fgas = b+a*xHI = (1-xHI)/zeta
    xHI(eq)=xHI1(eq); %changed
    currentIonization = 1-xHI;
    Mcrit=calcMcritV2Fast(z,zMax,N,zs,currentIonization,fin, ionizationHistory) ;
    ionizationHistory=[];


else
    Mcrit = zer;
end


end


