classdef SIM21Gets
    methods(Static)
        function Hz = getHz(z)
            % Hubble constant
            % in km/s/Mpc
            load(SIM21Utils.getMatrixPath('Planck_parameters'));
            Hz = H0*sqrt(Om*(1+z).^3+OLambda+8.5522e-05*(1+z).^4); 
        end
        
        
        function Kappa = getKappa(TK)
            TkLib = [1 2 4 6 8 10 15 20 25 30 40 50 60 70 80 90 100 200 300 500 700 1000 2000 3000 5000 7000 10000 10^5 10^6 10^7 10^8];
            KappaLib=[1.38e-13 1.43e-13 2.71e-13 6.6e-13 1.47e-12 2.88e-12 9.10e-12 1.78e-11 2.73e-11 3.67e-11 5.38e-11 6.86e-11 8.14e-11,...
                    9.25e-11 1.02e-10 1.11e-10 1.19e-10 1.75e-10 2.09e-10 2.56e-10 2.91e-10 3.31e-10 4.27e-10 4.97e-10 6.03e-10 6.87e-10,...
                    7.87e-10 1.9919e-9 4.9385e-9 1.2244e-8 3.0357e-8];
            Kappa=exp(interp1(log(TkLib),log(KappaLib),log(TK)));
        end
        
        
        function  [JA_alpha] = getLyA(TK, z0,ncube,fstar,flag,flagM,XeffTerm,Ispec,Reion,feedback,p,pop,FSfunc,photoheatingVersion)
            global delta_cube
            global pathname_Data1
            global ID
            
            load(SIM21Utils.getMatrixPath('Planck_parameters'));
            
            N=length(delta_cube);
            
            JA_alpha = zeros(N,N,N);
            Zrange =5:66; %9:50; 
            Ind1 = find(Zrange > z0,1,'first');  
            zint =  Zrange(Ind1-1:Ind1); 
            if (length(zint)==2)
                JA_interp = zeros(2,N,N,N);
                for kk=1:2      
                    if(zint(kk)>50)
                        JA_interp(kk,:,:,:) = (1e-20)*ones(N,N,N);%+JAX;
                    else
                        load(strcat(pathname_Data1,'Jalpha_',num2str(zint(kk)),ID,'.mat'));
                        %JA_interp(kk,:,:,:) = JA+(1+delta_cube.*LWgetDz( zint(kk))/LWgetDz(40)).*JAX*nb*(1+zint(kk))^3*c/(getHz(zint(kk))*4*pi*nu_a^2*hpl)/(3*10^24)^2; 
                        JA_interp(kk,:,:,:) = Jalpha;%+JAX;
                        Jalpha =[]; 
                    end
                    %JAX =[];
                end      
                JA_alpha = exp(squeeze((interp1(log(1+zint),log(JA_interp),log(1+z0)))));
                JA_interp=[];
            else
                load(strcat(pathname_Data1,'Jalpha_',num2str(zint),ID,'.mat'));
                JA_alpha =Jalpha;%X+JA;%
            end
            Jalpha=[];
            JAX=[];
        end
        
        
        function nH = getnH(z)
            % hydrogen proper number density
            load(SIM21Utils.getMatrixPath('Planck_parameters'));
            Y = 0.247; % Helium abundance by mass
            nH=(rhoc/mp)*(1-Y)*Ob*(1+z).^3; % 1/Mpc^3
        end
        
        
        function Sa = getSalpha(Tk,z)
            Sa=exp(-.803*(Tk).^(-2/3)*(1e-6*SIM21Gets.getTauGPn(z)).^(1/3));
        end
        
        
        function tGP = getTauGPn(z) 
            load(SIM21Utils.getMatrixPath('Planck_parameters'));
            gamma = 50*10^6; % 1/s
            nH = SIM21Gets.getnH(z);% 1/Mpc^3
            tGP=(3*nH*(lambda_a)^3*gamma)./(2*SIM21Gets.getHz(z)*3.241e-20);
        end
        
        
        function Ts = getTs(TK,z0,ncube,fstar,flag,flagM,XeffTerm,Ispec,Reion,feedback,p,pop,FSfunc,photoheatingVersion)
            % calculate spin temperature
            Tcmb = SIM21Gets.getTcmb(z0); 
            xA = SIM21Gets.getXA(TK,z0,ncube,fstar,flag,flagM,XeffTerm,Ispec,Reion,feedback,p,pop,FSfunc,photoheatingVersion);
            xC = SIM21Gets.getXc(TK,z0); 
            Tse = 0.4;%[K]
            
            Ts = (1+xC+xA.*(1+Tse./TK).^(-1))./(Tcmb^(-1)+xC.*TK.^(-1)+xA.*TK.^(-1).*(1+Tse./TK).^(-1));
        end
        
        
        function [x_talpha]=getXA(TK,z0,ncube,fstar,flag,flagM,XeffTerm,Ispec,Reion,feedback,p,pop,FSfunc,photoheatingVersion) 
            %  z0=16;
            %  ncube=0;
            %  fstar=0.05; 
            %  flag=0;
            %  flagM=16.5;
            %  XeffTerm=1;
            %  Ispec=1;
            %  Reion=0.075;
            %  feedback=0;
            %  p=0;
            load(SIM21Utils.getMatrixPath('Planck_parameters'));
            
            S_alpha2=SIM21Gets.getSalpha(TK,z0);
            JA_alpha = SIM21Gets.getLyA(TK,z0,ncube,fstar,flag,flagM,XeffTerm,Ispec,Reion,feedback,p,pop,FSfunc,photoheatingVersion);
            x_talpha=real(1.7e11*(1+z0).^(-1)*(S_alpha2.*JA_alpha));
        end
        
        
        function Xc = getXc(TK,z0)
            % collisional coupling coefficient
            global delta_cube
            load(SIM21Utils.getMatrixPath('Planck_parameters'));
            
            nH = (ones(size(delta_cube))+max(-1,min(1,delta_cube))).*SIM21Gets.getnH(z0)/(3.1e24)^3;% 1/cm^3
            kappa = SIM21Gets.getKappa(TK);%cm^3/s
            Tcmb = SIM21Gets.getTcmb(z0);%K
            Tstar = 0.068; %K
            Xc = (nH.*kappa/A10)*(Tstar/Tcmb);
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        function f = getf(z,A,B,Lion)   
            global delta_cube
            load(SIM21Utils.getMatrixPath('Planck_parameters'));
            
            dtdz = -3.09e19./(SIM21Gets.getHz(z).*(1+z)); % in seconds
            fHe = 1-0.921;
            fHI = 0.921;
            nb0 = 2.1e-7*(1+fHe/fHI);
            alphaA=4.2e-13*(exp(B)/10^4).^(-0.76);
            f = (1+z).*dtdz*Lion./exp(A)-(1+z)^4.*dtdz.*fHI*nb0.*alphaA.*((exp(A)-1).^2./exp(A)).*(1+delta_cube*LWgetDz(z)/D40);
        end
        
        
        function g = getg(z,A,B,Lion,eps)
            fHe=6/(6+76);
            kb = (1.38e-16/1.6e-12); % eV/K
            dtdz =-3.09e19./(SIM21Gets.getHz(z).*(1+z)); % in seconds
            g = 2-SIM21Gets.getf(z,A,B,Lion)+(2/3/kb)*(1+z).*dtdz*eps./(exp(B).*exp(A))+((exp(A)-1)./(exp(A)+fHe)).*8.55e-13.*(2.725*(1+z)-exp(B)).*(1+z).^5./(exp(B));
        end
        
        
        function Rc = getHorizon(zcenter)
            load(SIM21Utils.getMatrixPath('Planck_parameters'));
            F = @(x)c./(H0*sqrt(Om*(1+x).^3+OLambda));
            Rc = quadl(F,zcenter,1e10);
        end
        
        
        function zs = getRtoz(R,z)
            % output: z+dz

            H0=67.04;
            c=3e5;
            omm=.3169;
            k = 2*c/(H0*sqrt(omm));
            Rcrit = k./sqrt(1+z);
            
            if (R<SIM21Gets.getHorizon(z))
                zs = (k^2*z + 2*k*R*sqrt(1+z) - R.^2.*(1+z))./(k^2 - 2*k*R*sqrt(1+z) + R.^2*(1+z));
            else 
                zs = 1e15;
            end
        end
        
        
        function zmax = getzmax(zc,n)
            zmax=((1+zc)*((1-(n+1)^-2)/(1-n^-2)))-1;
        end
        
        
        function t = getztot(z)
            t = ((1+z)/7).^(-3/2)*.95e9;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        function Tcmb = getTcmb(z)
            Tcmb = 2.725*(1+z); 
        end
    end
end