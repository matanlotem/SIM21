function  [JA_alpha] = getLyA(TK, z0,ncube,fstar,flag,flagM,XeffTerm,Ispec,Reion,feedback,p,pop,FSfunc,photoheatingVersion)
    load Planck_parameters
    
    global delta_cube
    global pathname_Data1
    global ID
    
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