% 17 march 2014

% with Aviad
%power spectrum with velocity term added
% Reconstruction

function [Pk,Pk_iso,Pk_X,Pk_del,Pmu4,Pmu2,Pmu0,Kout,nk] = getPkRec(cube,del,Lx,ep)
%    N = length(cube);
%    dx = Lx/N;
%    pix = dx;
%    Vpix = dx^3;
%    
%    kx = zeros(size(cube));
%    ky = zeros(size(cube));
%    kz = zeros(size(cube));
%    
%    for Ind = 0:N-1
%        kx(Ind+1,:,:) = -(Ind>(N-1)/2)*(2*pi*(N-Ind))/(N*pix)*ones(1,N,N)+(Ind<=(N-1)/2)*(2*pi*Ind/(N*pix))*ones(1,N,N);
%        ky(:,Ind+1,:) = -(Ind>(N-1)/2)*(2*pi*(N-Ind))/(N*pix)*ones(N,1,N)+(Ind<=(N-1)/2)*(2*pi*Ind/(N*pix))*ones(N,1,N); 
%        kz(:,:,Ind+1) = -(Ind>(N-1)/2)*(2*pi*(N-Ind))/(N*pix)*ones(N,N,1)+(Ind<=(N-1)/2)*(2*pi*Ind/(N*pix))*ones(N,N,1);
%    end
%    
%    kk = (kx.^2+ky.^2+kz.^2).^0.5;
%    
%    %--- adding velocity term (20.01.2014) --------
%    mu = kz./kk;
%    Ind = kk==0;
%    size(Ind)
%    mu(Ind) = 0;
%    
%    k0=2*pi*((2*ep+1)/(2*ep))/Lx;
%    powK=0:100;
%    KN = k0*((2*ep+1)/(2*ep-1)).^powK;
%    
%    Kout = (2*pi/Lx)*((2*ep+1)/(2*ep-1)).^(1:100);
    %[kk,KN,Vpix,mu,Kout] = getPkVars(cube,Lx,ep);
    [kk,Vpix,mu] = SIM21Pk.getPkVars(cube,Lx);
    [KN,Kout] = SIM21Pk.getKout(Lx,ep);
    
    kdel = fftn(del);
    kcube = fftn(cube);
    kcube_del = kdel.*mu.^2;
    kcube_new = kcube+kcube_del;
    
    % OLD CODE
    %Pk = zeros(1,length(KN)-1); 
    %Pk_iso = zeros(1,length(KN)-1); 
    %Pk_X = zeros(1,length(KN)-1); 
    %Pk_del = zeros(1,length(KN)-1); 
    %nk = zeros(1,length(KN)-1); 
    
    %Pmu4= zeros(1,length(KN)-1);
    %Pmu2= zeros(1,length(KN)-1);
    %Pmu0= zeros(1,length(KN)-1);
    % END OLD CODE
    
    
    [NNN, bink] = histc(kk,KN);
    
    %MATAN NEW CODE
    lKN1=length(KN)-1;
    NMbink = repmat(bink,[1,1,1,lKN1]) == permute(repmat((1:lKN1)',[1,size(bink)]),[2,3,4,1]);
    Nnk = squeeze(sum(sum(sum(NMbink))))';


    Mkcn = repmat((kcube_new.*conj(kcube_new)),[1,1,1,lKN1]).*NMbink;
    Mkc = repmat((kcube.*conj(kcube)),[1,1,1,lKN1]).*NMbink;
    Mkcd = repmat((kcube_del.*conj(kcube_del)),[1,1,1,lKN1]).*NMbink;
    Mkckcd = repmat((kcube.*conj(kcube_del)),[1,1,1,lKN1]).*NMbink + repmat((kcube_del.*conj(kcube)),[1,1,1,lKN1]).*NMbink;
    Mkd = repmat((kdel.*conj(kdel)),[1,1,1,lKN1]).*NMbink;
    Mkckd = repmat((kcube.*conj(kdel)),[1,1,1,lKN1]).*NMbink + repmat((kdel.*conj(kcube)),[1,1,1,lKN1]).*NMbink;
    
    function N = calcN(M)
        N = 1./(Lx*Lx*Lx)*(Vpix^2)*squeeze(sum(sum(sum(M))))'./Nnk;
        N(Nnk==0) = 0;
    end
    NPk = calcN(Mkcn);
    NPk_iso = calcN(Mkc);
    NPk_X = calcN(Mkckcd);
    NPk_del = calcN(Mkcd);

    NPmu4 = calcN(Mkd);
    NPmu2 = calcN(Mkckd);
    NPmu0 = NPk_iso;
    
    nk=Nnk;
    Pk=NPk;
    Pk_iso=NPk_iso;
    Pk_X=NPk_X;
    Pk_del=NPk_del;
    Pmu4=NPmu4;
    Pmu2=NPmu2;
    Pmu0=NPmu0;
    
    % OLD CODE
    %for jk=1:length(KN)-1
    %    Mbink = bink==jk; 
    %    nk(jk) = sum(sum(sum(Mbink)));
    %    
    %    if(nk(jk)==0)
    %        Pk(jk) = 0;
    %        Pk_iso(jk) = 0;
    %        Pk_X(jk) = 0;
    %        Pk_del(jk) = 0;
    %    else
    %        Pk(jk) = 1./(Lx*Lx*Lx)*(Vpix)^2*sum(sum(sum((kcube_new).*conj(kcube_new).*Mbink)))./nk(jk);
    %        Pk_iso(jk) = 1./(Lx*Lx*Lx)*(Vpix)^2*sum(sum(sum((kcube).*conj(kcube).*Mbink)))./nk(jk);
    %        Pk_X(jk) = 1./(Lx*Lx*Lx)*(Vpix)^2*sum(sum(sum((kcube).*conj(kcube_del).*Mbink+(kcube_del).*conj(kcube).*Mbink)))./nk(jk);
    %        Pk_del(jk) = 1./(Lx*Lx*Lx)*(Vpix)^2*sum(sum(sum((kcube_del).*conj(kcube_del).*Mbink)))./nk(jk);
    %        
    %        Pmu4(jk) = 1./(Lx*Lx*Lx)*(Vpix)^2*sum(sum(sum((kdel).*conj(kdel).*Mbink)))./nk(jk); 
    %        Pmu2(jk) = 1./(Lx*Lx*Lx)*(Vpix)^2*sum(sum(sum((kcube).*conj(kdel).*Mbink+(kdel).*conj(kcube).*Mbink)))./nk(jk);
    %        Pmu0(jk) = 1./(Lx*Lx*Lx)*(Vpix)^2*sum(sum(sum((kcube).*conj(kcube).*Mbink)))./nk(jk);
    %    end
    %end
    % END OLD CODE
end