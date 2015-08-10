function [kk,KN,Vpix,mu,Kout] = getPkVars(cube,Lx,ep)
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
    
    k0=2*pi*((2*ep+1)/(2*ep))/Lx;
    powK=0:100;
    KN = k0*((2*ep+1)/(2*ep-1)).^powK;
    
    Kout = (2*pi/Lx)*((2*ep+1)/(2*ep-1)).^(1:100);
end