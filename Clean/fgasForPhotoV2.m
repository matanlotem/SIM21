function fgas = fgasForPhotoV2(flag,flagM,JLW21,z,fstarM,fstarA,FSfunc,currentIonization,Mcrit)
global vbc_cube
global delta_cube

if (z>15)
    delmax = 1.0;
    Nres = 100;
    deltas = linspace(-delmax,delmax,Nres);
    vbc = linspace(0,3.8,Nres);
    NresM = 45;
    if (z<66)
    Mass = logspace(log10(1e5),log10(1e11),NresM);
    else
    Mass = logspace(log10(1e5),log10(1e8),NresM);
    end

else
     
    delmax = 1.0;
    Nres = 100;
    deltas = linspace(-delmax,2*delmax,3*Nres/2);
    vbc = linspace(0,3.8,Nres);
    NresM = 45;
   Mass = logspace(log10(1e5),log10(1e11),NresM);
end

delta_c=1.686-1.686*(1+z)/2000;
dz = 0.0001;
g = LWgetDz(z)/LWgetDz(40);
gp = LWgetDz(z+dz)/LWgetDz(40);

        


    
McubeGM = LWgetMcubeVc(flag, z, JLW21,flagM);   

if(z>65)
    gridA=0;gridM=0;
    load(strcat('/scratch300/matanlotem/DataBackgrounds_withPlanck/gridM_',num2str(z),'.mat'));
    load(strcat('/scratch300/matanlotem/DataBackgrounds_withPlanck/gridA_',num2str(z),'.mat'));
else
     
    if FSfunc==1
    gridA=0;gridM=0;
    load(strcat('/scratch300/matanlotem/DataBackgrounds_withPlanck/gridM10_',num2str(z),'.mat'));
    load(strcat('/scratch300/matanlotem/DataBackgrounds_withPlanck/gridA10_',num2str(z),'.mat'));  
     end
      if FSfunc==2
    gridA=0;gridM=0;
    load(strcat('/scratch300/matanlotem/DataBackgrounds_withPlanck/gridM10Sharp_',num2str(z),'.mat'));
    load(strcat('/scratch300/matanlotem/DataBackgrounds_withPlanck/gridA10Sharp_',num2str(z),'.mat'));  
     end  
end


fgasA = exp(interp3(vbc,deltas,Mass,log(gridA),(flag+1e-10)*vbc_cube,...
    max(min(delta_cube*g,0.99*delta_c*ones(size(delta_cube))),-0.99*ones(size(delta_cube))),McubeGM,'linear'));

fgasM = exp(interp3(vbc,deltas,Mass,log(gridM),(flag+1e-10)*vbc_cube,...
    max(min(delta_cube*g,0.99*delta_c*ones(size(delta_cube))),-0.99*ones(size(delta_cube))),McubeGM,'linear'));

    fgas = fstarM*fgasM+fstarA*fgasA;
    
    
    
    McubeGM = max(McubeGM,Mcrit);

    fgasA = exp(interp3(vbc,deltas,Mass,log(gridA),(flag+1e-10)*vbc_cube,...
        max(min(delta_cube*g,0.99*delta_c*ones(size(delta_cube))),-0.99*ones(size(delta_cube))),McubeGM,'linear'));
    gridA=0;
    fgasM = exp(interp3(vbc,deltas,Mass,log(gridM),(flag+1e-10)*vbc_cube,...
        max(min(delta_cube*g,0.99*delta_c*ones(size(delta_cube))),-0.99*ones(size(delta_cube))),McubeGM,'linear'));
    gridM=0;
    fgasNew = fstarM*fgasM+fstarA*fgasA;

    fgas = fgas.*(1-currentIonization) + fgasNew.*currentIonization;

end



