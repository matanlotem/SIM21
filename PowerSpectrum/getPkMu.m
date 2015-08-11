%
function [Pk,Pkm,Kout,Muout,nk,nkm] = getPkMu(cube,Lx,ep,Nmu)
kcube = fftn(cube);
%kcube = fftshift(kcube);
N = length(kcube);
dx = Lx/N;
pix = dx;
Vpix = dx^3;

k0=2*pi*((2*ep+1)/(2*ep))/Lx;
powK=0:100;
KN = k0*((2*ep+1)/(2*ep-1)).^powK;

Kout = (2*pi/Lx)*((2*ep+1)/(2*ep-1)).^(1:100);

Muvec = linspace(-1,1,Nmu+1);
%del = 2/Nmu;
%Muout = linspace(-1+del,1-del,Nmu);
Muout = (Muvec(1:end-1)+Muvec(2:end))/2;
Pk = zeros(1,length(KN)-1); 
nk = zeros(1,length(KN)-1); 

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
Ind = find(kk==0);
mu(Ind) = 0;


[NNN, bink] = histc(kk,KN);


        for jk=1:length(KN)-1
            
              Mbink = bink==jk; 
              nk(jk) = sum(sum(sum(Mbink)));
              
              if(nk(jk)==0)
                 Pk(jk) = 0;
              else
               Pk(jk) = 1./(Lx*Lx*Lx)*(Vpix)^2*sum(sum(sum((kcube).*conj(kcube).*Mbink)))./nk(jk);
              end
           
        end

Pkm = zeros(Nmu,length(KN)-1); 
nkm = zeros(Nmu,length(KN)-1); 
[MMM, mink] = histc(mu,Muvec);
% the spectrum is symmetric
        for jk=1:length(KN)-1
            
              Mbink = bink==jk; 
              for jm=1:Nmu
                  Mmink = (mink==jm); 
                  nkm(jm,jk) = sum(sum(sum(Mbink.*Mmink)));

                  if(nkm(jm,jk)==0)
                     Pkm(jm,jk) = 0;
                  else
                   Pkm(jm,jk) = 1./(Lx*Lx*Lx)*(Vpix)^2*sum(sum(sum((kcube).*conj(kcube).*Mbink.*Mmink)))./nkm(jm,jk);
                  end
              end
           
        end
end