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

function [temp_adj,pres_adj,psal_adj,time_adj,datamode] = profiles_adjusted(dir_name)

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

temp_adj = NaN(1100,i); pres_adj = NaN(1100,i); psal_adj = NaN(1100,i); 
datamode = char(NaN(i,1)); time_adj = NaN(i,1);  

i = 1;

for n = 1:totalfiles
    
    disp(n);
    
    fn = num2str(float_num(n));
    
    filename = strcat(dir_name,'Profiles/',fn,'_prof.nc');
    [t,p,s,~,~,~,jd,dm] = read_prof_adj(filename);
    
    if float_num(n) == 2900171
        disp('2900171');
%         tqc = [tqc(1:28,:); tqc(30:end,:)]; tqc = tqc(:,1:50);
%         sqc = [sqc(1:28,:); sqc(30:end,:)]; sqc = sqc(:,1:50);
%         pqc = [pqc(1:28,:); pqc(30:end,:)]; pqc = pqc(:,1:50);
%         posqc = [posqc(1:28); posqc(30:end)];
%         juldqc = [juldqc(1:28); juldqc(30:end)];
%         proftqc = [proftqc(1:28); proftqc(30:end)];
%         profsqc = [profsqc(1:28); profsqc(30:end)];
%         profpqc = [profpqc(1:28); profpqc(30:end)];
    end
    
    [mm,nn] = size(p);
    
    [ii,jj] = meshgrid((i:i-1+nn),(1:mm));
    
    ind = sub2ind(size(psal_adj),jj,ii);
   
    if ~isempty(s)
    psal_adj(ind) = s;
    end 
    
    if ~isempty(t)
    temp_adj(ind) = t;
    end 
    
    if ~isempty(p)
    pres_adj(ind) = p;
    end 
    
    datamode(i:i-1+nn,1) = dm;
    time_adj(i:i-1+nn,1) = jd;
    i = i+nn;
    
    clear t p s terr perr serr dm
end

% save /home/alison/Matlab/Mat' files'/Global' Mapping 2'/North' Atlantic'/profiles_adj.mat
save /home/alison/Matlab/Mat' files'/Global' Mapping 2'/South' Atlantic'/profiles_adj_SA.mat
% save /home/alison/Matlab/Mat' files'/Global' Mapping 2'/Indian/profiles_adj.mat
% save /home/alison/Matlab/Mat' files'/Global' Mapping 2'/North' Pacific'/profiles_adj.mat
% save /home/alison/Matlab/Mat' files'/Global' Mapping 2'/South' Pacific'/profiles_adj.mat


return




