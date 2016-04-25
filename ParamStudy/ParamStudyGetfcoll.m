function fcollTOT=ParamStudyGetfcoll(runCase,z)
    pathname_DataBackgrounds = '/scratch300/matanlotem/DataBackgrounds_withPlanck/';
    vbc_cube=importdata(strcat(pathname_DataBackgrounds,'my',num2str(runCase.ncube),'_v.dat'));

    vbc = mean(mean(mean(vbc_cube)));

    if runCase.feedback
    	LWm = mean(mean(mean(runCase.getData('JLW',z,1))));
    else
    	LWm = 0;
    end

    flag = 1;

    %LWgetMcubeVc
    h = 0.6704;
    Oc = 0.12038/h^2;
    Ob = 0.022032/h^2;
    Om = Ob+Oc;
    Omz = Om*(1+z).^3./(Om.*(1+z).^3+1-Om);
    d = Omz-1;
    Dc = 18*pi^2+82*d-39*d^2;

    a=4.015;
    zrec = 1020;
    Vc = (runCase.vc^2 + (a*flag*vbc*(1+z)/(1+zrec)*0.000097*3e5).^2).^0.5;

    Mc = 10^8*(Vc/23.4).^3*((1+z)/10).^(-3/2)*(Om*Dc/Omz/18/pi^2).^-0.5/h;
    M = Mc.*(1+ 6.96*(LWm*4*pi).^0.47);

    M = min(M,1e11*ones(size(M)));
    M = max(M,10^5*ones(size(M)));


    if runCase.fsfunc == 1
        gridA=0;gridM=0;    
        gridAMQ=0;gridMMQ=0;
        load([pathname_DataBackgrounds,'gridM10_',num2str(SIM21Utils.roundDataTypeZ('xHI',z)),'.mat']);
        load([pathname_DataBackgrounds,'gridA10_',num2str(SIM21Utils.roundDataTypeZ('xHI',z)),'.mat']);  
    elseif runCase.fsfunc == 2
        gridA=0;gridM=0;
        gridAMQ=0;gridMMQ=0;
        load([pathname_DataBackgrounds,'gridM10Sharp_',num2str(SIM21Utils.roundDataTypeZ('xHI',z)),'.mat']);
        load([pathname_DataBackgrounds,'gridA10Sharp_',num2str(SIM21Utils.roundDataTypeZ('xHI',z)),'.mat']); 
    end

    delmax = 1.0;
    Nres = 100;
    if (z>15)
        deltas = linspace(-delmax,delmax,Nres);
    else
        deltas = linspace(-delmax,2*delmax,3*Nres/2);
    end
    delta=deltas(51);
    vbc_vals = linspace(0,3.8,Nres);
    NresM = 45;
    if (z<66)
        Mass = logspace(log10(1e5),log10(1e11),NresM);
    else
        Mass = logspace(log10(1e5),log10(1e8),NresM);
    end


    %%%%%%%%%%%%%
    zz = SIM21Utils.getDataType('xHI').z;
    zz = zz(zz<=50);
    
    zIND=find(zz<z,1,'last');
    xHImat(zIND:length(zz)) = mean(mean(mean(runCase.getData('xHI',zz(zIND:end)),2),3),4);

    ZxHI = interp1(zz(zIND:end),xHImat(zIND:end),z);
    q = find(xHImat>(1-(1-ZxHI)/2));
    ZhalfxHI = interp1(xHImat(q-1:q),zz(q-1:q),1-(1-ZxHI)/2);


    M0 = 2.8e9; 
    a = 0.17;
    b = -2.1;
    c = 2;
    d = 2.5;
    J21 = 0.5;  

    Mcrit = M0.*(J21.^a).*(((1+z)./10).^b).*((1-((1 + z)./(1 + ZhalfxHI)).^c).^d);
    Mmin=max(M,Mcrit);
    %%%%%%%%%%%%%%%

    fgasA = exp(interp3(vbc_vals,deltas,Mass,log(gridA),(flag+1e-10)*vbc,delta,M,'linear'));
    fgasM = exp(interp3(vbc_vals,deltas,Mass,log(gridM),(flag+1e-10)*vbc,delta,M,'linear'));
    fcoll=fgasM+fgasA;

    fgasA = exp(interp3(vbc_vals,deltas,Mass,log(gridA),(flag+1e-10)*vbc,delta,Mmin,'linear'));
    fgasM = exp(interp3(vbc_vals,deltas,Mass,log(gridM),(flag+1e-10)*vbc,delta,Mmin,'linear'));
    fcollTOT=(1-ZxHI)*(fgasM+fgasA)+ZxHI*fcoll;
end