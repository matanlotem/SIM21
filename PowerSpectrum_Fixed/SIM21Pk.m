classdef SIM21Pk
    methods(Static)
        function [kk,Vpix,mu] = getPkVars(cube,Lx)
            N = length(cube);
            dx = Lx/N;
            pix = dx;
            Vpix = dx^3;
            
            kx = zeros(size(cube));
            ky = zeros(size(cube));
            kz = zeros(size(cube));
            
            for Ind = 0:N-1
                kx(Ind+1,:,:) = -(Ind>(N-1)/2)*(2*pi*(N-Ind))/(N*pix)*ones(1,N,N)+(Ind<=(N-1)/2)*(2*pi*Ind/(N*pix))*ones(1,N,N);
                ky(:,Ind+1,:) = -(Ind>(N-1)/2)*(2*pi*(N-Ind))/(N*pix)*ones(N,1,N)+(Ind<=(N-1)/2)*(2*pi*Ind/(N*pix))*ones(N,1,N); 
                kz(:,:,Ind+1) = -(Ind>(N-1)/2)*(2*pi*(N-Ind))/(N*pix)*ones(N,N,1)+(Ind<=(N-1)/2)*(2*pi*Ind/(N*pix))*ones(N,N,1);
            end
            
            kk = (kx.^2+ky.^2+kz.^2).^0.5;
            
            %--- adding velocity term (20.01.2014) --------
            mu = kz./kk;
            %Ind = find(kk==0); %Ind = kk==0; %%% I think this is a bug - in getPkRec line 30
            Ind = kk==0;
            mu(Ind) = 0;
        end
        
        
        function [KN,Kout] = getKout(Lx,ep)
            k0=2*pi*((2*ep+1)/(2*ep))/Lx;
            powK=0:100;
            KN = k0*((2*ep+1)/(2*ep-1)).^powK;
            
            Kout = (2*pi/Lx)*((2*ep+1)/(2*ep-1)).^(1:100);
        end
        
        
        function [delta_out] = getDV(del,Lx)
            [kk,Vpix,mu] = SIM21Pk.getPkVars(del,Lx);
            kdel = fftn(del);
            kcube_del = kdel.*mu.^2;
            delta_out = (ifftn(kcube_del));
        end
        
        
        function [Pk,Pk_iso,Pk_X,Pk_del,Pmu4,Pmu2,Pmu0,Kout,nk] = getPkRec(cube,del,Lx,ep)
            [kk,Vpix,mu] = SIM21Pk.getPkVars(cube,Lx);
            [KN,Kout] = SIM21Pk.getKout(Lx,ep);
            
            kdel = fftn(del);
            kcube = fftn(cube);
            kcube_del = kdel.*mu.^2;
            kcube_new = kcube+kcube_del;
            
            [NNN, bink] = histc(kk,KN);
            
            lKN1=length(KN)-1;
            NMbink = repmat(bink,[1,1,1,lKN1]) == permute(repmat((1:lKN1)',[1,size(bink)]),[2,3,4,1]);
            nk = squeeze(sum(sum(sum(NMbink))))';
            
            function out = calcThis(tmp_cube)
                M = repmat(tmp_cube,[1,1,1,lKN1]).*NMbink;
                out = 1./(Lx*Lx*Lx)*(Vpix^2)*squeeze(sum(sum(sum(M))))'./nk;
                out(nk==0) = 0;
            end
            
            Pk = calcThis(kcube_new.*conj(kcube_new));
            Pk_iso = calcThis(kcube.*conj(kcube));
            Pk_X = calcThis((kcube.*conj(kcube_del))+(kcube_del.*conj(kcube)));
            Pk_del = calcThis(kcube_del.*conj(kcube_del));
            
            Pmu4 = calcThis(kdel.*conj(kdel));
            Pmu2 = calcThis((kcube.*conj(kdel)) + (kdel.*conj(kcube)));
            Pmu0 = Pk_iso;
        end
        
        
        %function [Pk,Pkm,Kout,Muout,nk,nkm] = getPkMu(cube,Lx,ep,Nmu)
        function [Pkm,Kout,Muout,nkm] = getPkMu(cube,Lx,ep,Nmu)
            [kk,Vpix,mu] = SIM21Pk.getPkVars(cube,Lx);
            [KN,Kout] = SIM21Pk.getKout(Lx,ep);
            
            Muvec = linspace(-1,1,Nmu+1);
            Muout = (Muvec(1:end-1)+Muvec(2:end))/2;
            
            kcube = fftn(cube);
            
            [NNN, bink] = histc(kk,KN);
            [MMM, mink] = histc(mu,Muvec);
            
            lKN1=length(KN)-1;
            NMbink = repmat(bink,[1,1,1,lKN1]) == permute(repmat((1:lKN1)',[1,size(bink)]),[2,3,4,1]);
            NMmink = repmat(mink,[1,1,1,Nmu]) == permute(repmat((1:Nmu)',[1,size(mink)]),[2,3,4,1]);
            
            %%%%%!!!!!!!!!!!! REQUIRES A LOT OF MEMORY !!!!!!!!!!!!!
            NMbm = repmat(NMbink,[1,1,1,1,Nmu]).*permute(repmat(NMmink,[1,1,1,1,lKN1]),[1,2,3,5,4]);
            
            %%%%%% WHY? This is calculated elsewhere...
            %nk = squeeze(sum(sum(sum(NMbink))))';
            %Pk = 1./(Lx*Lx*Lx)*(Vpix^2)*squeeze(sum(sum(sum(repmat((kcube.*conj(kcube)),[1,1,1,lKN1]).*NMbink))))'./nk;
            %Pk(nk==0) = 0;
            %%%%%%
            NMbink = [];
            NMmink = [];
            
            nkm = squeeze(sum(sum(sum((NMbm)))))';
            NMbmk = squeeze(sum(sum(sum(repmat(kcube.*conj(kcube),[1,1,1,lKN1,Nmu]).*NMbm))))';
            
            Pkm = 1./(Lx*Lx*Lx)*(Vpix^2)*NMbmk./nkm;
            Pkm(nkm==0) = 0;
        end
    end
end