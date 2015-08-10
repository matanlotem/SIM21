function Mcrit=calcMcritV2Fast(z,zMax,N,zs,currentIonization,fin, ionizationHistory)  

previousZIndex=find((round(zs*10)/10)==(round(z*10)/10))-1;
halfCurrentIonization = currentIonization.*0.5;
currentIonization=[];

MinHalf=min(min(min(halfCurrentIonization)));

%[ionizationHistory,fin]=getIonHis(directory,filePrefix,N,z,zs,MinHalf);






       %% Build a NxNxNxD matrix of halfCurrentIonization matrixes (D such)
    comparator = repmat(halfCurrentIonization,[1,1,1,fin]);
    

    %% Build a boolean matrix where 1 is more than half the current ionization
    ionizationLarger = logical(ionizationHistory(:,:,:,1:fin) >= comparator);
    
    
    comparator=[];
    
    for i=2:fin
        ionizationLarger(:,:,:,i)=or(ionizationLarger(:,:,:,i),ionizationLarger(:,:,:,i-1));
    end

    %% Find the first index where ionization was larger than half the current
    firstIndexes = fin - sum(ionizationLarger,4) + 1;
    
   


    A = (sum(ionizationLarger,4)==0);
    firstIndexes(A)=fin;
    
    
    zIN=zs(firstIndexes-fin+previousZIndex);
    
  


    ind = bsxfun(@plus, bsxfun(@plus,...
        (1:N).', N*(0:N-1)), N*N*permute(0:N-1, [3 1 2])) + N*N*N*(firstIndexes-1);
    zINIonization = ionizationHistory(ind);
  

    prevFirstIndexes = firstIndexes-1;
    prevFirstIndexesIsZero = logical(prevFirstIndexes==0);
    prevFirstIndexes(prevFirstIndexesIsZero) = 1;

    ind = bsxfun(@plus, bsxfun(@plus,...
        (1:N).', N*(0:N-1)), N*N*permute(0:N-1, [3 1 2])) + N*N*N*((prevFirstIndexes)-1);
    zINPrevIonization = ionizationHistory(ind);

    ionizationHistory=[];
    
    zINPrev = zs(prevFirstIndexes-fin+previousZIndex);

    
    zIN = 10.^(((log10(zINIonization)-log10(halfCurrentIonization)).*log10(zINPrev) +...
        (log10(halfCurrentIonization)-log10(zINPrevIonization)).*log10(zIN))./(log10(zINIonization)-log10(zINPrevIonization)));
    %zIN(doNotInterpolateHere) = lastzIN(doNotInterpolateHere);

    %%
    %% Calculate M crit
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Only update M crit from 0 where the latest history is at least half ionized
    %% (ionizationLarger(:,:,:,previousZIndex) is not defined yet)
    if z+2<=zMax
        zIN_ion = ionizationLarger(:,:,:,fin);
        ionizationLarger=[];
        
    else
        zIN_ion = true([N,N,N]);
    end


    %% Best fit parameters
    M0 = 2.8e9; % In solar masses
    a = 0.17;
    b = -2.1;
    c = 2;
    d = 2.5;
    J21 = 0.5;  % Factor of two over 0.25, the value for J in HII regions
                %  for all z (at all redshifts).
    Mcrit = zeros(N,N,N);
    Mcrit(zIN_ion) = M0.*(J21.^a).*(((1+z)./10).^b).*((1-((1 + z)./(1 + zIN(zIN_ion))).^c).^d);
    


end