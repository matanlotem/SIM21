%% 
%% Calculates the critical mass for the photoheating effect (version 1)
%% + directory - Where the xHI files are stored.
%% + filePrefix - A.K.A. '_xHI_'.
%% + zMax - Maximal value of z in the simulation for which the function
%%          will be called upon.
%% + N - Mass matrix dimension size.
%% + z - The redshift to calculate for.
%%
%% To run the tests set DEBUG=1 (define globally)
function Mcrit = calculateMcritV1changes(directory, filePrefix, zMax, N, z)

    global DEBUG;
    if isempty(DEBUG)
        DEBUG = 0;
    end
    global gridInterpKey;

    %persistent zIN; % Carefull using persistent this here
                    % , its designed for a single simulation calling this function
    %persistent zIN_ion;
    %if isempty(zIN) || isempty(zIN_ion)
    zIN = zeros(N,N,N);
    zIN_ion = false(N,N,N);
    %end
    %persistent nextZ;
    %if isempty(nextZ)
        nextZ = zMax;
    %end
    zs=[60:-1:16 15:-0.1:5];
    if z>=15
        previousZIndex = zMax-(z+1)+1; % Number of previous z
    end
    if z<15
        previousZIndex = round((zMax-(15+1)+1)+(15-z)*10); % Number of previous z
    end

    if z<zMax

        %%
        %% Calculate zIN and zIN_ion
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %% Load all saved files up to current z (Not efficient, we
        %% can define xHI globally. But its not terribly slow)
        IonizationHistory=zeros(previousZIndex,N,N,N);
        for i=1:previousZIndex

            %% Load xHI from file if needed
            load([directory filePrefix num2str(zs(i)) gridInterpKey]);
            IonizationHistory(i,:,:,:) = 1-xHI;

            %% ionized: 1 means is the cell is fully ionized, 0 otherwise
            ionized = logical(xHI==0);

            %% Find the new ionizations:
            %% (all the ones ionized in either zIN_ion or ionized but not both) and
            %% (in ionized)
            zIN_ion_new = xor(zIN_ion,ionized) & ionized;

            %% Mark all entries that are ionized
            zIN_ion = (ionized);

            %% Set newly ionized z as i
            zIN(zIN_ion_new) = zs(i);
            NewZ=zeros(N,N,N);
            if i>2
        
                PrevIon=zeros(N,N,N);
                Prev2Ion=zeros(N,N,N);
                PrevIon=squeeze(IonizationHistory(i-1,:,:,:));
                Prev2Ion=squeeze(IonizationHistory(i-2,:,:,:));

                %NewZ=squeeze(10.^interp1(log10(interp_ion+1e-60),log10([zs(i-1),zs(i-2)]),log10(ones(N,N,N)),'linear','extrap'));
                NewZ=10.^(((log10(Prev2Ion)-log10(ones(N,N,N)))*log10(zs(i-1))+(log10(ones(N,N,N))-log10(PrevIon))*log10(zs(i-2)))./(log10(Prev2Ion)-log10(PrevIon)));
                qq=or(NewZ>zs(i-1)-1e-5,NewZ<z+1e-5);
                NewZ(qq)=0;

            end

            zIN(zIN_ion_new)=max(zIN(zIN_ion_new),NewZ(zIN_ion_new));
               
        end
        if z<58
            PrevIon=zeros(N,N,N);
            Prev2Ion=zeros(N,N,N);
            PrevIon=squeeze(IonizationHistory(previousZIndex,:,:,:));
            Prev2Ion=squeeze(IonizationHistory(previousZIndex-1,:,:,:));
            %NewZ=squeeze(10.^interp1(log10(interp_ion+1e-60),log10([zs(i-1),zs(i-2)]),log10(ones(N,N,N)),'linear','extrap'));
            
            NewZ=10.^(((log10(Prev2Ion)-log10(ones(N,N,N)))*log10(zs(previousZIndex))+(log10(ones(N,N,N))-log10(PrevIon))*log10(zs(previousZIndex-1)))./(log10(Prev2Ion)-log10(PrevIon)));
            
            just_ionized = and(~(ionized),and(logical(z<NewZ),logical(NewZ<zs(previousZIndex))));
            
            zIN(just_ionized) = max(zIN(just_ionized), NewZ(just_ionized));
            zIN_ion = or(ionized,just_ionized);
        end
        
   
        %% Start from this z next time
        nextZ = min(z,zMax);

        %%
        %% Calculate Mcrit
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %% Best fit parameters
        M0 = 2.8e9; % In solar masses
        a = 0.17;
        b = -2.1;
        c = 2;
        d = 2.5;
        J21 = 0.5;  % Factor of two over 0.25, the value for J in HII regions
                    %  for all z (at all redshifts).
        %% Set M crit only for the ones ionized (the rest are left as 0)
        Mcrit = zeros(N,N,N);
        Mcrit(zIN_ion) = M0.*(J21.^a).*(((1+z)./10).^b).*((1-((1 + z)./(1 + zIN(zIN_ion))).^c).^d);

    else
        Mcrit = zeros(N,N,N);
    end
end