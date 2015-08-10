% 17 march 2014

% with Aviad
%power spectrum with velocity term added
% Reconstruction

function [delta_out] = getDV(del,Lx)
    kdel = fftn(del);
    
    N = length(kdel);
    dx = Lx/N;
    pix = dx;
    Vpix = dx^3;
    
    kx = zeros(size(del));
    ky = zeros(size(del));
    kz = zeros(size(del));
    
    for Ind = 0:N-1
        kx(Ind+1,:,:) = -(Ind>(N-1)/2)*(2*pi*(N-Ind))/(N*pix)*ones(1,N,N)+(Ind<=(N-1)/2)*(2*pi*Ind/(N*pix))*ones(1,N,N);
        ky(:,Ind+1,:) = -(Ind>(N-1)/2)*(2*pi*(N-Ind))/(N*pix)*ones(N,1,N)+(Ind<=(N-1)/2)*(2*pi*Ind/(N*pix))*ones(N,1,N); 
        kz(:,:,Ind+1) = -(Ind>(N-1)/2)*(2*pi*(N-Ind))/(N*pix)*ones(N,N,1)+(Ind<=(N-1)/2)*(2*pi*Ind/(N*pix))*ones(N,N,1);
    end
    
    kk = (kx.^2+ky.^2+kz.^2).^0.5;
    
    %--- adding velocity term (20.01.2014) --------
    mu = kz./kk; 
    Ind = kk==0;
    mu(Ind) = 0;
    
    kcube_del = kdel.*mu.^2;
    
    delta_out = (ifftn(kcube_del));
end