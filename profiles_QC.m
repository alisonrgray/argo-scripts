%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%
% user-set variables - examples:%
%
% dir_name = '/home/alison/Research/Data/Argo Profiles/Float/North Atlantic/';
% dir_name = '/home/alison/Data/Argo Profiles/Float/South Atlantic/';
% dir_name = '/home/alison/Research/Data/Argo Profiles/Float/Indian/';
% dir_name = '/home/alison/Research/Data/Argo Profiles/Float/North Pacific/';
% dir_name = '/home/alison/Research/Data/Argo Profiles/Float/South Pacific/';
%%%%%%%%%%%%%%%%%%%%%%

function [temp_qc,pres_qc,psal_qc,position_qc,juld_qc,prof_temp_qc,prof_psal_qc,prof_pres_qc] = profiles_QC(dir_name)

cd(strcat(dir_name,'Profiles'));
direc = dir;
floats = {};
[floats{1:length(direc),1}] = deal(direc.name);
floats = strtok(floats,'_');
floats = floats(3:end,:);

float_num = zeros(length(floats),1);

for i = 1:length(floats);
    float_num(i) = str2double(char(floats(i)));
end

clear direc floats i

cd('/home/alison/Matlab/M files/Global Mapping')
totalfiles = length(float_num);

i = 1;
for n = 1:totalfiles
    disp(n);
    
    fn = num2str(float_num(n));
    
    filename = strcat(dir_name,'Profiles/',fn,'_prof.nc');
    [~,p] = read_prof(filename);
    
    if float_num(n) == 2900171
        p = [p(:,1:28) p(:,30:end)]; p = p(1:50,:);
    end
    
    [~,nn] = size(p);
    
    i = i+nn;
    
    clear p nn mm filename
    
end

position_qc = char(NaN(i,1));  juld_qc = char(NaN(i,1));
temp_qc = char(NaN(i,1100)); pres_qc = char(NaN(i,1100)); psal_qc = char(NaN(i,1100));
prof_temp_qc = char(NaN(i,1));  prof_psal_qc = char(NaN(i,1));
prof_pres_qc = char(NaN(i,1));
temp_qc_adj = char(NaN(i,1100)); pres_qc_adj = char(NaN(i,1100)); psal_qc_adj = char(NaN(i,1100));

i = 1;

for n = 1:totalfiles
    
    disp(n);
    
    fn = num2str(float_num(n));
    
    filename = strcat(dir_name,'Profiles/',fn,'_prof.nc');
    [tqc,pqc,sqc,posqc,juldqc,proftqc,profsqc,profpqc,tqc_adj,pqc_adj,sqc_adj] = read_prof_qc(filename);
    tqc = tqc'; pqc = pqc'; sqc = sqc';
    tqc_adj = tqc_adj'; pqc_adj = pqc_adj'; sqc_adj = sqc_adj';
    
    
    if float_num(n) == 2900171
        tqc = [tqc(1:28,:); tqc(30:end,:)]; tqc = tqc(:,1:50);
        sqc = [sqc(1:28,:); sqc(30:end,:)]; sqc = sqc(:,1:50);
        pqc = [pqc(1:28,:); pqc(30:end,:)]; pqc = pqc(:,1:50);
        posqc = [posqc(1:28); posqc(30:end)];
        juldqc = [juldqc(1:28); juldqc(30:end)];
        proftqc = [proftqc(1:28); proftqc(30:end)];
        profsqc = [profsqc(1:28); profsqc(30:end)];
        profpqc = [profpqc(1:28); profpqc(30:end)];
    end
    
    [mm,nn] = size(pqc);
    
    position_qc(i:i-1+mm) = posqc;
    juld_qc(i:i-1+mm) = juldqc;
    
    if ~isempty(sqc)
        psal_qc(i:i-1+mm,1:nn) = sqc;
        psal_qc_adj(i:i-1+mm,1:nn) = sqc_adj;
        prof_psal_qc(i:i-1+mm) = profsqc;
    end
    
    if ~isempty(tqc)
        temp_qc(i:i-1+mm,1:nn) = tqc;
        temp_qc_adj(i:i-1+mm,1:nn) = tqc_adj;
        prof_temp_qc(i:i-1+mm) = proftqc;
    end
    
    if ~isempty(pqc)
        pres_qc(i:i-1+mm,1:nn) = pqc;
        pres_qc_adj(i:i-1+mm,1:nn) = pqc_adj;
        prof_pres_qc(i:i-1+mm) = profpqc;
    end
    
    i = i+mm;
    
    clear tqc pqc sqc posqc juldqc proftqc profsqc profpqc pqc_adj tqc_adj sqc_adj
end

clear dir_name filename float_num fn i mm n nn totalfiles

% save /home/alison/Matlab/Mat' files'/Global' Mapping'/North' Atlantic'/profiles_QC.mat
save /home/alison/Matlab/Mat' files'/Global' Mapping 2'/South' Atlantic'/profiles_QC_SA.mat
% save /home/alison/Matlab/Mat' files'/Global' Mapping'/Indian/profiles_QC.mat
% save /home/alison/Matlab/Mat' files'/Global' Mapping'/North' Pacific'/profiles_QC.mat
% save /home/alison/Matlab/Mat' files'/Global' Mapping'/South' Pacific'/profiles_QC.mat


return




