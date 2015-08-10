pathname_Output = '/scratch300/matanlotem/ParamStudy/';
TFiles = dir([pathname_Output,'P*']);
for i = 1:length(TFiles)
    TFile = TFiles(i);
    flag = isequal(importdata([pathname_Output,TFile.name]),importdata([pathname_Output,'N_',TFile.name]));
    disp([num2str(flag),' ',TFile.name]);
end