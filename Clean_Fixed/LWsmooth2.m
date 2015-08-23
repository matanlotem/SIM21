% uploads window from LWmakeWindow.m
function sdensityf= LWsmooth2(kdensityf,R)
    global pathname_DataBackgrounds
    %load(strcat('/home/anastasia/LWproject/my_work/My_LW/WR_',num2str(R),'.mat')); 
    %load(strcat('/scratch300/Anastasia_F/DataBackgrounds/WR_',num2str(R),'.mat')); 
    if (length(kdensityf)==256)
        load(strcat(pathname_DataBackgrounds,'WR256_',num2str(R),'.mat')); 
    else
        load(strcat(pathname_DataBackgrounds,'WR128_',num2str(R),'.mat'));
    end
        
    sdensityf = real(ifftn(windowk.*squeeze(kdensityf)));

end

