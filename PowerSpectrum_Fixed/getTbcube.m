function [Tb,T21] = getTbcube(zii,ncube,fstar,flag,flagM,XeffTerm,Ispec,Reion,feedback,p,pop,FSfunc,photoheatingVersion)
    global pathname_Data2
    global ID
    global delta_cube
    load Planck_parameters
    
    del = real(max(-0.9,min(1,delta_cube*LWgetDz(zii)/LWgetDz(40))));
    
    if(ncube<100)
        Lpix = 3;
    else
        Lpix = 6;
    end
    
    N=length(delta_cube);
    ep=10;
    Lx=Lpix*N;
    C21 = 27*sqrt(0.15/Om/h^2)*(Ob*h^2/0.023);
    z21cm=5:50;
    zion = [5:0.1:15, 16:1:50];
    
    Ind1 = find(z21cm>zii,1,'first');  
    zint = z21cm(Ind1-1:Ind1);  
    F_interp = zeros(2,N,N,N);
    
    for iz=1:2
        load(strcat(pathname_Data2,'TK_',num2str(zint(iz)),ID,'.mat'));
        Ts = SIM21Gets.getTs(TK,zint(iz),ncube,fstar,flag,flagM,XeffTerm,Ispec,Reion,feedback,p,pop,FSfunc,photoheatingVersion);
        %Ts = TK;
        Ts(Ts==0)=1e-20;   
        deltaTerm = (ones(size(delta_cube))+max(-0.9,min(1,delta_cube*LWgetDz(zint(iz))/LWgetDz(40))));
        
        Tb = C21*sqrt((1+zint(iz))/10).*deltaTerm.*(Ts-2.725*(1+zint(iz)))./Ts;
        F_interp(iz,:,:,:) = Tb;
        Tb =[];        
    end
    
    Tb = squeeze((interp1(log(1+zint),F_interp,log(1+zii)))); 
    F_interp=[];
    
    load(strcat(pathname_Data2,'xHI_',num2str(zii),ID,'.mat'));
    load(strcat(pathname_Data2,'Neut_',num2str(zii),ID,'.mat'));
    
    Tb =  Tb.*xHI.*Neut;
    T21 = mean(mean(mean(Tb)));
   
end