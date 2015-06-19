%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%
% user-set variables - examples:%
% 
% dir_name = '/home/alison/Data/Argo Profiles/Float/North Pacific/';
% dir_name = '/home/alison/Research/Data/Argo Profiles/Float/South Pacific/';
% dir_name = '/home/alison/Data/Argo Profiles/Float/North Atlantic/';
% dir_name = '/home/alison/Data/Argo Profiles/Float/South Atlantic/';
% dir_name = '/home/alison/Data/Argo Profiles/Float/Indian/';
%%%%%%%%%%%%%%%%%%%%%%

function [longitude,latitude,time,temp,pres,psal,fnum,type] = profiles(dir_name)

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

cd('/home/alison/Matlab/M files/Global Mapping 2')
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

longitude = NaN(i,1); latitude = NaN(i,1);  time = NaN(i,1); 
temp = NaN(1150,i); pres = NaN(1150,i); psal = NaN(1150,i); 
fnum = NaN(i,1); type = NaN(i,1);

i = 1;
for n = 1:totalfiles
    
    disp(n);
    
    fn = num2str(float_num(n));
    
    filename = strcat(dir_name,'Profiles/',fn,'_prof.nc');
    [t,p,s,ln,lt,jd,ty] = read_prof(filename);
   
    if float_num(n) == 2900171
        t = [t(:,1:28) t(:,30:end)]; t = t(1:50,:);
        s = [s(:,1:28) s(:,30:end)]; s = s(1:50,:);
        p = [p(:,1:28) p(:,30:end)]; p = p(1:50,:);
        lt = [lt(1:28); lt(30:end)];
        ln = [ln(1:28); ln(30:end)];
        jd = [jd(1:28); jd(30:end)];
        ty = [ty(1:28); ty(30:end)];
    end
    
    [mm,nn] = size(p);
    
    [ii,jj] = meshgrid((i:i-1+nn),(1:mm));
    
    ind = sub2ind(size(psal),jj,ii);
    
    longitude(i:i-1+nn) = ln;
    latitude(i:i-1+nn) = lt;
    time(i:i-1+nn) = jd;
   
    if ~isempty(s)
    psal(ind) = s;
    end 
    
    if ~isempty(t)
    temp(ind) = t;
    end 
    
    if ~isempty(p)
    pres(ind) = p;
    end 
    
    fnum(i:i-1+nn) = float_num(n)*ones(nn,1);
    
    if length(ty) == length(jd)
    type(i:i-1+nn) = ty;
    else if length(ty) == 3
            type(i:i-1+nn) = (str2num(strcat(num2str(ty(1)),num2str(ty(2)),num2str(ty(3)))))*ones(nn,1); %#ok<ST2NM>
        end
    end
    i = i+nn;
    
    clear t p s ln lt jd ty mm nn ii jj ind

end

clear ans dir_name filename fn i n totalfiles

% save /home/alison/Matlab/Mat' files'/Global' Mapping'/North' Atlantic'/profiles.mat
% save /home/alison/Matlab/Mat' files'/Global' Mapping'/South' Atlantic'/profiles.mat
% save /home/alison/Matlab/Mat' files'/Global' Mapping'/Indian/profiles.mat
% save /home/alison/Matlab/Mat' files'/Global' Mapping'/North' Pacific'/profiles.mat
% save /home/alison/Matlab/Mat' files'/Global' Mapping'/South' Pacific'/profiles.mat

% save /home/alison/Matlab/Mat' files'/Global' Mapping 2'/South' Atlantic'/profiles_SA.mat
% save /home/alison/Matlab/Mat' files'/Global' Mapping 2'/North' Atlantic'/profiles_NA.mat
% save /home/alison/Matlab/Mat' files'/Global' Mapping 2'/North' Pacific'/profiles_NP.mat
save /home/alison/Matlab/Mat' files'/Global' Mapping 2'/Indian/profiles_IND.mat -v7.3

return


    
    
    
    
    
    