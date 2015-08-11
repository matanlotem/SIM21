% 17 march 2014

% with Aviad
%power spectrum with velocity term added
% Reconstruction

function [Pk,Pk_iso,Pk_X,Pk_del,Pmu4,Pmu2,Pmu0,Kout,nk] = getPkRec(cube,del,Lx,ep)

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

kcube = fftn(cube);
kcube_del = kdel.*mu.^2;
kcube_new = kcube+kcube_del;
%---------------------------------


k0=2*pi*((2*ep+1)/(2*ep))/Lx;
powK=0:100;
KN = k0*((2*ep+1)/(2*ep-1)).^powK;
Kout = (2*pi/Lx)*((2*ep+1)/(2*ep-1)).^(1:100);
Pk = zeros(1,length(KN)-1); 
Pk_iso = zeros(1,length(KN)-1); 
Pk_X = zeros(1,length(KN)-1); 
Pk_del = zeros(1,length(KN)-1); 
nk = zeros(1,length(KN)-1); 

Pmu4= zeros(1,length(KN)-1);
Pmu2= zeros(1,length(KN)-1);
Pmu0= zeros(1,length(KN)-1);


[NNN, bink] = histc(kk,KN);

        for jk=1:length(KN)-1
            
              Mbink = bink==jk; 
              nk(jk) = sum(sum(sum(Mbink)));
              
              if(nk(jk)==0)
                 Pk(jk) = 0;
                 Pk_iso(jk) = 0;
                 Pk_X(jk) = 0;
                 Pk_del(jk) = 0;
              else
               Pk(jk) = 1./(Lx*Lx*Lx)*(Vpix)^2*sum(sum(sum((kcube_new).*conj(kcube_new).*Mbink)))./nk(jk);
               Pk_iso(jk) = 1./(Lx*Lx*Lx)*(Vpix)^2*sum(sum(sum((kcube).*conj(kcube).*Mbink)))./nk(jk);
               Pk_X(jk) = 1./(Lx*Lx*Lx)*(Vpix)^2*sum(sum(sum((kcube).*conj(kcube_del).*Mbink+(kcube_del).*conj(kcube).*Mbink)))./nk(jk);
               Pk_del(jk) = 1./(Lx*Lx*Lx)*(Vpix)^2*sum(sum(sum((kcube_del).*conj(kcube_del).*Mbink)))./nk(jk);
               
               Pmu4(jk) = 1./(Lx*Lx*Lx)*(Vpix)^2*sum(sum(sum((kdel).*conj(kdel).*Mbink)))./nk(jk); 
               Pmu2(jk) = 1./(Lx*Lx*Lx)*(Vpix)^2*sum(sum(sum((kcube).*conj(kdel).*Mbink+(kdel).*conj(kcube).*Mbink)))./nk(jk);
               Pmu0(jk) = 1./(Lx*Lx*Lx)*(Vpix)^2*sum(sum(sum((kcube).*conj(kcube).*Mbink)))./nk(jk);
              end
           
        end
   
    
    
end