%% North Pacific
% 2005
clear
pack

load /home/alison/Matlab/Mat' files'/Global' Mapping 2'/North' Pacific'/profiles_NP.mat
load /home/alison/Matlab/Mat' files'/Global' Mapping 2'/North' Pacific'/profiles_QC_NP.mat
t1 = datenum(2005,01,01);
t2 = datenum(2006,01,01);
ind = find(time >= t1 & time <= t2);

temp = temp(:,ind); pres = pres(:,ind); psal = psal(:,ind);
latitude = latitude(ind); longitude = longitude(ind); time = time(ind);
fnum = fnum(ind); type = type(ind);

temp_qc = temp_qc(ind,:); pres_qc = pres_qc(ind,:); psal_qc = psal_qc(ind,:);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
temp_qc_adj = temp_qc_adj(ind,:); pres_qc_adj = pres_qc_adj(ind,:); psal_qc_adj = psal_qc_adj(ind,:);

clear t1 t2 dir_name ind float_num

save('/home/alison/Matlab/Mat files/Global Mapping 2/North Pacific/prof_2005.mat');

clear; pack;
load /home/alison/Matlab/Mat' files'/Global' Mapping 2'/North' Pacific'/profiles_adj_NP.mat
t1 = datenum(2005,01,01);
t2 = datenum(2006,01,01);
ind = find(time_adj >= t1 & time_adj <= t2);

temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind); time_adj = time_adj(ind);
datamode = datamode(ind);

clear t1 t2 dir_name ind float_num

save('/home/alison/Matlab/Mat files/Global Mapping 2/North Pacific/prof_adj_2005.mat');

clear; pack;
load /home/alison/Matlab/Mat' files'/Global' Mapping 2'/North' Pacific'/profiles_err_NP.mat
load /home/alison/Matlab/Mat' files'/Global' Mapping 2'/North' Pacific'/profiles_adj_NP.mat time_adj
t1 = datenum(2005,01,01);
t2 = datenum(2006,01,01);
ind = find(time_adj >= t1 & time_adj <= t2);

temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind); time_adj = time_adj(ind);

clear t1 t2 dir_name ind float_num

save('/home/alison/Matlab/Mat files/Global Mapping 2/North Pacific/prof_err_2005.mat');

% 2006	
clear
pack

load /home/alison/Matlab/Mat' files'/Global' Mapping 2'/North' Pacific'/profiles_NP.mat
load /home/alison/Matlab/Mat' files'/Global' Mapping 2'/North' Pacific'/profiles_QC_NP.mat
t1 = datenum(2006,01,01);
t2 = datenum(2007,01,01);
ind = find(time >= t1 & time <= t2);

temp = temp(:,ind); pres = pres(:,ind); psal = psal(:,ind);
latitude = latitude(ind); longitude = longitude(ind); time = time(ind);
fnum = fnum(ind); type = type(ind);

temp_qc = temp_qc(ind,:); pres_qc = pres_qc(ind,:); psal_qc = psal_qc(ind,:);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
temp_qc_adj = temp_qc_adj(ind,:); pres_qc_adj = pres_qc_adj(ind,:); psal_qc_adj = psal_qc_adj(ind,:);

clear t1 t2 dir_name ind float_num

save('/home/alison/Matlab/Mat files/Global Mapping 2/North Pacific/prof_2006.mat');

clear; pack;
load /home/alison/Matlab/Mat' files'/Global' Mapping 2'/North' Pacific'/profiles_adj_NP.mat
t1 = datenum(2006,01,01);
t2 = datenum(2007,01,01);
ind = find(time_adj >= t1 & time_adj <= t2);

temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind); time_adj = time_adj(ind);
datamode = datamode(ind);

clear t1 t2 dir_name ind float_num

save('/home/alison/Matlab/Mat files/Global Mapping 2/North Pacific/prof_adj_2006.mat');

clear; pack;
load /home/alison/Matlab/Mat' files'/Global' Mapping 2'/North' Pacific'/profiles_err_NP.mat
load /home/alison/Matlab/Mat' files'/Global' Mapping 2'/North' Pacific'/profiles_adj_NP.mat time_adj
t1 = datenum(2006,01,01);
t2 = datenum(2007,01,01);
ind = find(time_adj >= t1 & time_adj <= t2);

temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind); time_adj = time_adj(ind);

clear t1 t2 dir_name ind float_num

save('/home/alison/Matlab/Mat files/Global Mapping 2/North Pacific/prof_err_2006.mat');

% 2007
clear
pack

load /home/alison/Matlab/Mat' files'/Global' Mapping 2'/North' Pacific'/profiles_NP.mat
load /home/alison/Matlab/Mat' files'/Global' Mapping 2'/North' Pacific'/profiles_QC_NP.mat
t1 = datenum(2007,01,01);
t2 = datenum(2008,01,01);
ind = find(time >= t1 & time <= t2);

temp = temp(:,ind); pres = pres(:,ind); psal = psal(:,ind);
latitude = latitude(ind); longitude = longitude(ind); time = time(ind);
fnum = fnum(ind); type = type(ind);

temp_qc = temp_qc(ind,:); pres_qc = pres_qc(ind,:); psal_qc = psal_qc(ind,:);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
temp_qc_adj = temp_qc_adj(ind,:); pres_qc_adj = pres_qc_adj(ind,:); psal_qc_adj = psal_qc_adj(ind,:);

clear t1 t2 dir_name ind float_num

save('/home/alison/Matlab/Mat files/Global Mapping 2/North Pacific/prof_2007.mat');

clear; pack;
load /home/alison/Matlab/Mat' files'/Global' Mapping 2'/North' Pacific'/profiles_adj_NP.mat
t1 = datenum(2007,01,01);
t2 = datenum(2008,01,01);
ind = find(time_adj >= t1 & time_adj <= t2);

temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind); time_adj = time_adj(ind);
datamode = datamode(ind);

clear t1 t2 dir_name ind float_num

save('/home/alison/Matlab/Mat files/Global Mapping 2/North Pacific/prof_adj_2007.mat');

clear; pack;
load /home/alison/Matlab/Mat' files'/Global' Mapping 2'/North' Pacific'/profiles_err_NP.mat
load /home/alison/Matlab/Mat' files'/Global' Mapping 2'/North' Pacific'/profiles_adj_NP.mat time_adj
t1 = datenum(2007,01,01);
t2 = datenum(2008,01,01);
ind = find(time_adj >= t1 & time_adj <= t2);

temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind); time_adj = time_adj(ind);

clear t1 t2 dir_name ind float_num

save('/home/alison/Matlab/Mat files/Global Mapping 2/North Pacific/prof_err_2007.mat');

% 2008
clear
pack

load /home/alison/Matlab/Mat' files'/Global' Mapping 2'/North' Pacific'/profiles_NP.mat
load /home/alison/Matlab/Mat' files'/Global' Mapping 2'/North' Pacific'/profiles_QC_NP.mat
t1 = datenum(2008,01,01);
t2 = datenum(2009,01,01);
ind = find(time >= t1 & time <= t2);

temp = temp(:,ind); pres = pres(:,ind); psal = psal(:,ind);
latitude = latitude(ind); longitude = longitude(ind); time = time(ind);
fnum = fnum(ind); type = type(ind);

temp_qc = temp_qc(ind,:); pres_qc = pres_qc(ind,:); psal_qc = psal_qc(ind,:);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
temp_qc_adj = temp_qc_adj(ind,:); pres_qc_adj = pres_qc_adj(ind,:); psal_qc_adj = psal_qc_adj(ind,:);

clear t1 t2 dir_name ind float_num

save('/home/alison/Matlab/Mat files/Global Mapping 2/North Pacific/prof_2008.mat');

clear; pack;
load /home/alison/Matlab/Mat' files'/Global' Mapping 2'/North' Pacific'/profiles_adj_NP.mat
t1 = datenum(2008,01,01);
t2 = datenum(2009,01,01);
ind = find(time_adj >= t1 & time_adj <= t2);

temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind); time_adj = time_adj(ind);
datamode = datamode(ind);

clear t1 t2 dir_name ind float_num

save('/home/alison/Matlab/Mat files/Global Mapping 2/North Pacific/prof_adj_2008.mat');

clear; pack;
load /home/alison/Matlab/Mat' files'/Global' Mapping 2'/North' Pacific'/profiles_err_NP.mat
load /home/alison/Matlab/Mat' files'/Global' Mapping 2'/North' Pacific'/profiles_adj_NP.mat time_adj
t1 = datenum(2008,01,01);
t2 = datenum(2009,01,01);
ind = find(time_adj >= t1 & time_adj <= t2);

temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind); time_adj = time_adj(ind);

clear t1 t2 dir_name ind float_num

save('/home/alison/Matlab/Mat files/Global Mapping 2/North Pacific/prof_err_2008.mat');

% 2009
clear
pack

load /home/alison/Matlab/Mat' files'/Global' Mapping 2'/North' Pacific'/profiles_NP.mat
load /home/alison/Matlab/Mat' files'/Global' Mapping 2'/North' Pacific'/profiles_QC_NP.mat
t1 = datenum(2009,01,01);
t2 = datenum(2010,01,01);
ind = find(time >= t1 & time <= t2);

temp = temp(:,ind); pres = pres(:,ind); psal = psal(:,ind);
latitude = latitude(ind); longitude = longitude(ind); time = time(ind);
fnum = fnum(ind); type = type(ind);

temp_qc = temp_qc(ind,:); pres_qc = pres_qc(ind,:); psal_qc = psal_qc(ind,:);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
temp_qc_adj = temp_qc_adj(ind,:); pres_qc_adj = pres_qc_adj(ind,:); psal_qc_adj = psal_qc_adj(ind,:);

clear t1 t2 dir_name ind float_num

save('/home/alison/Matlab/Mat files/Global Mapping 2/North Pacific/prof_2009.mat');

clear; pack;
load /home/alison/Matlab/Mat' files'/Global' Mapping 2'/North' Pacific'/profiles_adj_NP.mat
t1 = datenum(2009,01,01);
t2 = datenum(2010,01,01);
ind = find(time_adj >= t1 & time_adj <= t2);

temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind); time_adj = time_adj(ind);
datamode = datamode(ind);

clear t1 t2 dir_name ind float_num

save('/home/alison/Matlab/Mat files/Global Mapping 2/North Pacific/prof_adj_2009.mat');

clear; pack;
load /home/alison/Matlab/Mat' files'/Global' Mapping 2'/North' Pacific'/profiles_err_NP.mat
load /home/alison/Matlab/Mat' files'/Global' Mapping 2'/North' Pacific'/profiles_adj_NP.mat time_adj
t1 = datenum(2009,01,01);
t2 = datenum(2010,01,01);
ind = find(time_adj >= t1 & time_adj <= t2);

temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind); time_adj = time_adj(ind);

clear t1 t2 dir_name ind float_num

save('/home/alison/Matlab/Mat files/Global Mapping 2/North Pacific/prof_err_2009.mat');

% 2010
clear
pack

load /home/alison/Matlab/Mat' files'/Global' Mapping 2'/North' Pacific'/profiles_NP.mat
load /home/alison/Matlab/Mat' files'/Global' Mapping 2'/North' Pacific'/profiles_QC_NP.mat
t1 = datenum(2010,01,01);
t2 = datenum(2011,01,01);
ind = find(time >= t1 & time <= t2);

temp = temp(:,ind); pres = pres(:,ind); psal = psal(:,ind);
latitude = latitude(ind); longitude = longitude(ind); time = time(ind);
fnum = fnum(ind); type = type(ind);

temp_qc = temp_qc(ind,:); pres_qc = pres_qc(ind,:); psal_qc = psal_qc(ind,:);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
temp_qc_adj = temp_qc_adj(ind,:); pres_qc_adj = pres_qc_adj(ind,:); psal_qc_adj = psal_qc_adj(ind,:);

clear t1 t2 dir_name ind float_num

save('/home/alison/Matlab/Mat files/Global Mapping 2/North Pacific/prof_2010.mat');

clear; pack;
load /home/alison/Matlab/Mat' files'/Global' Mapping 2'/North' Pacific'/profiles_adj_NP.mat
t1 = datenum(2010,01,01);
t2 = datenum(2011,01,01);
ind = find(time_adj >= t1 & time_adj <= t2);

temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind); time_adj = time_adj(ind);
datamode = datamode(ind);

clear t1 t2 dir_name ind float_num

save('/home/alison/Matlab/Mat files/Global Mapping 2/North Pacific/prof_adj_2010.mat');

clear; pack;
load /home/alison/Matlab/Mat' files'/Global' Mapping 2'/North' Pacific'/profiles_err_NP.mat
load /home/alison/Matlab/Mat' files'/Global' Mapping 2'/North' Pacific'/profiles_adj_NP.mat time_adj
t1 = datenum(2010,01,01);
t2 = datenum(2011,01,01);
ind = find(time_adj >= t1 & time_adj <= t2);

temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind); time_adj = time_adj(ind);

clear t1 t2 dir_name ind float_num

save('/home/alison/Matlab/Mat files/Global Mapping 2/North Pacific/prof_err_2010.mat');




%% quality control

% 2005
clear; pack
load /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/prof_2005.mat;
load /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/prof_adj_2005.mat;
load /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/prof_err_2005.mat;

pres_qc = pres_qc';
psal_qc = psal_qc';
temp_qc = temp_qc';
pres_qc_adj = pres_qc_adj';
psal_qc_adj = psal_qc_adj';
temp_qc_adj = temp_qc_adj';

longitude(longitude < 0) = longitude(longitude < 0) + 360;

% plot locations of data
figure(1); clf
axesm('mercator','MapLatLimit',[-5 80],'MapLonLimit',[100 280])
plotm(latitude,longitude,'.');
geoshow('landareas.shp','FaceColor',[.9 .9 .9])
setm(gca, 'grid','on','meridianlabel','on','parallellabel','on','plinelocation',5,'mlinelocation',5,'plabellocation',5,'mlabellocation',5)
axis tight

figure(2); clf
plot(longitude,latitude,'.'); grid on;

% lat/lon range
ind = find(longitude <= 290 & longitude >= 100 & latitude <= 80 & latitude >= 0);
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

% remove gulf of mexico & atlantic profiles
ind = find(~(longitude >= 270 & latitude >= 16));
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

% remove sea of japan profiles
ind = find(~((longitude < 138 & latitude < 44 & latitude > 34.5) | (longitude < 141 & latitude < 44 & latitude > 37))); 
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

% plot locations of data
figure(1); clf
axesm('mercator','MapLatLimit',[-5 80],'MapLonLimit',[-100 50])
plotm(latitude,longitude,'.');
geoshow('landareas.shp','FaceColor',[.9 .9 .9])
setm(gca, 'grid','on','meridianlabel','on','parallellabel','on','plinelocation',10,'mlinelocation',10,'plabellocation',10,'mlabellocation',10)
axis tight

% look at maximum pressures
maxp = (max(abs(pres),[],1))';

% find & remove profiles with maxp < 800 db
indp = find(maxp < 800);
per800 = 100*length(indp)/length(maxp) %#ok<*NOPTS,*NASGU>

num_800 = [fnum(indp) time(indp)];

ind = find(maxp >= 800);
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

clear ans ind indp per800 t2 

save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_800
clear num_800 maxp

% delete all Nan profiles
clear ind; k = 1; 
for i = 1:length(fnum)
    t = temp(:,i); s = psal(:,i); p = pres(:,i);
    ss = sum(~isnan(s)); tt = sum(~isnan(t)); pp = sum(~isnan(p));
    
    if ~(tt == 0 || ss == 0 || pp == 0)
        ind(k) = i; k = k+1; %#ok<NAGROW>
    end
end
clear t s p ss tt pp k i
ind = ind';
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

indt = find(sum(~isnan(temp),2) == 0); %#ok<MXFND>
indp = find(sum(~isnan(pres),2) == 0); %#ok<MXFND>
inds = find(sum(~isnan(psal),2) == 0); %#ok<MXFND>
ind = max([min(indt) min(indp) min(inds)]);
temp = temp(1:ind,:); psal = psal(1:ind,:); pres = pres(1:ind,:);
temp_qc = temp_qc(1:ind,:); pres_qc = pres_qc(1:ind,:); psal_qc = psal_qc(1:ind,:);
temp_adj = temp_adj(1:ind,:); psal_adj = psal_adj(1:ind,:); pres_adj = pres_adj(1:ind,:);
temp_qc_adj = temp_qc_adj(1:ind,:); pres_qc_adj = pres_qc_adj(1:ind,:); psal_qc_adj = psal_qc_adj(1:ind,:);
temp_err = temp_err(1:ind,:); psal_err = psal_err(1:ind,:); pres_err = pres_err(1:ind,:);
clear indt indp inds ind

save /home/alison/Matlab/M'at files'/Global' Mapping 2'/North' Pacific'/tmp.mat

% graylist floats
cd('/home/alison/Matlab/M files/Global Mapping 2/');
[grayfloat,~,start_date,end_date,errtype] = graylist('/home/alison/Matlab/Mat files/Global Mapping 2/North Pacific/');


%% %%%% quality control flags %%%%

% examine flags on temperature data
[mm,nn] = size(temp);
tqc = NaN(size(temp));
for j = 1:nn
    for i = 1:mm
        if ~isnan(temp(i,j))
            if (strcmp(temp_qc_adj(i,j),'1') || strcmp(temp_qc_adj(i,j),'0') || strcmp(temp_qc_adj(i,j),'2') || strcmp(temp_qc_adj(i,j),' ') || double(temp_qc_adj(i,j))==0)
                tqc(i,j) = 1;
            else
                tqc(i,j) = 0;
            end
        end
    end
end

% examine flags on salinity data
[mm,nn] = size(psal);
sqc = NaN(size(psal));
for j = 1:nn
    for i = 1:mm
        if ~isnan(psal(i,j))
            if (strcmp(psal_qc_adj(i,j),'1') || strcmp(psal_qc_adj(i,j),'0') || strcmp(psal_qc_adj(i,j),'2') || strcmp(psal_qc_adj(i,j),' ') || double(psal_qc_adj(i,j))==0)
                sqc(i,j) = 1;
            else
                sqc(i,j) = 0;
            end
        end
    end
end

% examine flags on pressure data
[mm,nn] = size(pres);
pqc = NaN(size(pres));
for j = 1:nn
    for i = 1:mm
        if ~isnan(pres(i,j))
            if (strcmp(pres_qc_adj(i,j),'1') || strcmp(pres_qc_adj(i,j),'0') || strcmp(pres_qc_adj(i,j),'2') || strcmp(pres_qc_adj(i,j),' ') || double(pres_qc_adj(i,j))==0)
                pqc(i,j) = 1;
            else
                pqc(i,j) = 0;
            end
        end
    end
end

clear mm nn i j


% display number of profiles with one or more bad qc flag
num_bad_temp = length(find(sum(tqc == 0,1) > 0))
num_bad_pres = length(find(sum(pqc == 0,1) > 0))
num_bad_psal = length(find(sum(sqc == 0,1) > 0))
clear num_bad_*

% make 2nd variable set with data with qc ~= 0,1,2,' ' set to NaN
temp2 = temp_adj; temp2(tqc == 0) = NaN;
psal2 = psal_adj; psal2(sqc == 0) = NaN;
pres2 = pres_adj; pres2(pqc == 0) = NaN;

% need t, s, and p at every point - make NaN's consistent across all data
ind = find(isnan(pres2));
temp2(ind) = NaN;
psal2(ind) = NaN;

ind = find(isnan(psal2));
temp2(ind) = NaN;
pres2(ind) = NaN;

ind = find(isnan(temp2));
pres2(ind) = NaN;
psal2(ind) = NaN;

% delete profiles that are all NaNs
clear ind indb; k = 1; l = 1;
for i = 1:length(fnum)
    t = temp2(:,i); s = psal2(:,i); p = pres2(:,i);
    ss = sum(~isnan(s)); tt = sum(~isnan(t)); pp = sum(~isnan(p));
    
    if ~(tt == 0 || ss == 0 || pp == 0)
        ind(k) = i; k = k+1; %#ok<SAGROW>
    else
        indb(l) = i; l = l+1; %#ok<SAGROW>
    end
end
clear t s p ss tt pp k i l
ind = ind'; indb = indb';

num_badQC = [fnum(indb) time(indb)];
save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_badQC -APPEND
clear indb num_badQC

temp2 = temp2(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres2 = pres2(:,ind);
psal2 = psal2(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
pqc = pqc(:,ind); pres = pres(:,ind); psal = psal(:,ind); sqc = sqc(:,ind); temp = temp(:,ind); tqc = tqc(:,ind);
pres_adj = pres_adj(:,ind); pres_err = pres_err(:,ind); pres_qc = pres_qc(:,ind);
psal_adj = psal_adj(:,ind); psal_err = psal_err(:,ind); psal_qc = psal_qc(:,ind);
temp_adj = temp_adj(:,ind); temp_err = temp_err(:,ind); temp_qc = temp_qc(:,ind);
datamode = datamode(ind); time_adj = time_adj(ind);

clear ind

% find & remove profiles with maxp < 800 db
maxp = (max(abs(pres2),[],1))';
indp = find(maxp < 800);
per800 = 100*length(indp)/length(maxp) %#ok<*NOPTS,*NASGU>

num_800_2 = [fnum(indp) time(indp)];

ind = find(maxp >= 800);
temp2 = temp2(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres2 = pres2(:,ind);
psal2 = psal2(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
pqc = pqc(:,ind); pres = pres(:,ind); psal = psal(:,ind); sqc = sqc(:,ind); temp = temp(:,ind); tqc = tqc(:,ind);
pres_adj = pres_adj(:,ind); pres_err = pres_err(:,ind); pres_qc = pres_qc(:,ind);
psal_adj = psal_adj(:,ind); psal_err = psal_err(:,ind); psal_qc = psal_qc(:,ind);
temp_adj = temp_adj(:,ind); temp_err = temp_err(:,ind); temp_qc = temp_qc(:,ind);
datamode = datamode(ind); time_adj = time_adj(ind);
clear ind

clear ans ind indp per800 t2 

save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_800_2
clear num_800_2 maxp

% % compute percentage of good values in each profile
% [~,nn] = size(pres2);
% percentNaN = NaN(nn,3);
% for i = 1:nn
%     p = pqc(:,i); t = tqc(:,i); s = sqc(:,i);
%     np = sum(p == 0); nt = sum(t == 0); ns = sum(s == 0); 
%     npp = sum(~isnan(p)); ntt = sum(~isnan(t)); nss = sum(~isnan(s));
%     percentNaN(i,:) = 100*[np./npp nt/ntt ns/nss];
% end
% clear i p np lastp mm nn npp t nt ntt s ns nss
% 
% ind = find(percentNaN(:,1) > 25 | percentNaN(:,2) > 25 | percentNaN(:,3) > 25);
% for ii = 1:length(ind)
%     fnum(ind(ii))
%     figure(1); clf
%     plot(temp(:,ind(ii)),-pres(:,ind(ii)),'m.-'); hold on;
%     plot(temp2(:,ind(ii)),-pres2(:,ind(ii)),'r.-');
%     disp(temp_qc(:,ind(ii))');
%     disp(pres_qc(:,ind(ii))');
%     disp(psal_qc(:,ind(ii))');
%     figure(2); clf
%     plot(psal(:,ind(ii)),-pres(:,ind(ii)),'c.-'); hold on;
%     plot(psal2(:,ind(ii)),-pres2(:,ind(ii)),'b.-');
%     pause
% end
% clear ii ans
% 
% num_percentNaN = [fnum(ind) time(ind)];
% save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_percentNaN -APPEND
% clear num_percentNaN
% 
% temp2(:,ind) = NaN; pres2(:,ind) = NaN; psal2(:,ind) = NaN;

% find & delete short profiles
n = sum(~isnan(temp2),1)';
% for i = 1:10
%     disp(i); ind = find(n == i); length(ind)
%     figure(1); clf
%     plot(temp2(:,ind),-pres2(:,ind),'.-');
%     pause
% end

ind = find(n <= 10 & n > 0); length(ind)
num_shortprof = [fnum(ind) time(ind)];
save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_shortprof -APPEND
clear num_shortprof

temp2(:,ind) = NaN; pres2(:,ind) = NaN; psal2(:,ind) = NaN;
clear n i ind ans

% delete profiles that are all NaNs
clear ind indb; k = 1; l = 1;
for i = 1:length(fnum)
    t = temp2(:,i); s = psal2(:,i); p = pres2(:,i);
    ss = sum(~isnan(s)); tt = sum(~isnan(t)); pp = sum(~isnan(p));
    
    if ~(tt == 0 || ss == 0 || pp == 0)
        ind(k) = i; k = k+1; %#ok<SAGROW>
    else
        indb(l) = i; l = l+1; %#ok<SAGROW>
    end
end
clear t s p ss tt pp k i l
ind = ind'; indb = indb';

temp2 = temp2(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres2 = pres2(:,ind);
psal2 = psal2(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
pqc = pqc(:,ind); pres = pres(:,ind); psal = psal(:,ind); sqc = sqc(:,ind); temp = temp(:,ind); tqc = tqc(:,ind);
pres_adj = pres_adj(:,ind); pres_err = pres_err(:,ind); pres_qc = pres_qc(:,ind);
psal_adj = psal_adj(:,ind); psal_err = psal_err(:,ind); psal_qc = psal_qc(:,ind);
temp_adj = temp_adj(:,ind); temp_err = temp_err(:,ind); temp_qc = temp_qc(:,ind);
datamode = datamode(ind); time_adj = time_adj(ind);
clear ind indb 

save /home/alison/Matlab/M'at files'/Global' Mapping 2'/North' Pacific'/tmp.mat

%% remove greylist floats
k = 1;
badfloats = 0; type_err = cell(1,1); bad_start = 0; bad_end = 0;
for i = 1:length(grayfloat)
    if sum(grayfloat(i) == fnum) > 0
    if grayfloat(i) ~= badfloats
        badfloats(k,1) = grayfloat(i); %#ok<SAGROW>
        type_err(k,1) = errtype(i);
        bad_start(k,1) = start_date(i); %#ok<SAGROW>
        bad_end(k,1) = end_date(i); %#ok<SAGROW>
        k = k+1;
    end
    end
end
clear k i

% ind = find(~strcmp(type_err,'potential problem with bin assignment'));
% badfloats = badfloats(ind);
% type_err = type_err(ind);
% bad_start = bad_start(ind);
% bad_end = bad_end(ind);

badstart = cellstr(num2str(bad_start));
badend = cellstr(num2str(bad_end));

for i = 1:length(badstart)
    bs = char(badstart(i));
    bs = datenum([sscanf(bs(1:4),'%f'),sscanf(bs(5:6),'%f'),sscanf(bs(7:8),'%f')]);
    bad_st(i,1) = bs; %#ok<SAGROW>
    be = char(badend(i));
    if ~strcmp(be,'     NaN') && ~strcmp(be,'NaN')
        be = datenum([sscanf(be(1:4),'%f'),sscanf(be(5:6),'%f'),sscanf(be(7:8),'%f')]);
        bad_end(i) = be;
    else
        bad_end(i) = NaN;
    end
end

fnumb = []; timeb = [];
for j = 1:length(badfloats)
    if isnan(bad_end(j))
        ind = find(~((fnum == badfloats(j))&(time > bad_st(j))));
        indb = find((fnum == badfloats(j)) & (time > bad_st(j)));
    else
        ind = find(~(fnum == badfloats(j) & time > bad_st(j) & time < bad_end(j)));
        indb = find(fnum == badfloats(j) & time > bad_st(j) & time < bad_end(j));
    end
fnumb = [fnumb; fnum(indb)]; %#ok<AGROW>
timeb = [timeb; time(indb)]; %#ok<AGROW>

temp2 = temp2(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres2 = pres2(:,ind);
psal2 = psal2(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
pqc = pqc(:,ind); pres = pres(:,ind); psal = psal(:,ind); sqc = sqc(:,ind); temp = temp(:,ind); tqc = tqc(:,ind);
pres_adj = pres_adj(:,ind); pres_err = pres_err(:,ind); pres_qc = pres_qc(:,ind);
psal_adj = psal_adj(:,ind); psal_err = psal_err(:,ind); psal_qc = psal_qc(:,ind);
temp_adj = temp_adj(:,ind); temp_err = temp_err(:,ind); temp_qc = temp_qc(:,ind);
datamode = datamode(ind); time_adj = time_adj(ind);

end

num_gray = [fnumb timeb];
save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_gray -APPEND

clear num_gray fnumb timeb indb ind
clear bad* end_date j ind param start_date be bs i 
clear ans errtype grayfloat ii percentNaN type_err

% graylist floats
cd('/home/alison/Matlab/M files/Global Mapping 2/');
[grayfloat,param,start_date,end_date,errtype] = graylist('/home/alison/Matlab/Mat files/Global Mapping 2/North Pacific/');

save /home/alison/Matlab/M'at files'/Global' Mapping 2'/North' Pacific'/tmp.mat

%% look for out-of-range pressures
figure(1); clf; imagesc(pres); colorbar;
figure(2); clf; imagesc(pres2); colorbar

outofrangenum = [];

ind = find(pres > 2200); length(ind)
[i,j] = find(pres2 > 2200); length(i)
tempnum = fnum(j);
presnum = 0; k = 1;
for ii = 1:length(tempnum)
if tempnum(ii) ~= presnum
    presnum(k,1) = tempnum(ii); %#ok<SAGROW>
    k = k+1;
end
end
outofrangenum = [outofrangenum; presnum];

% for ii = 1:length(i)
%     figure(1); clf
%     plot(temp2(i(ii),:),-pres2(i(ii),:),'r.-'); hold on
%     plot(psal2(i(ii),:),-pres2(i(ii),:),'b.-');
%     pres2(i(ii),1:30)
%     pause
% end

ind = find(pres < 0); length(ind)
[~,j] = find(pres2 < 0); length(j)
tempnum = fnum(j);
presnum = 0; k = 1;
for ii = 1:length(tempnum)
if tempnum(ii) ~= presnum
    presnum(k,1) = tempnum(ii); %#ok<SAGROW>
    k = k+1;
end
end
outofrangenum = [outofrangenum; presnum];

% for ii = 1:length(i)
%     figure(1); clf
%     plot(temp(i(ii),:),-pres(i(ii),:),'r.-'); hold on
%     plot(psal(i(ii),:),-pres(i(ii),:),'b.-');
%     pres2(i(ii),1:30)
%     pause
% end

% look for out-of-range temperatures
figure(1); clf; imagesc(temp); colorbar;
figure(2); clf; imagesc(temp2); colorbar

ind = find(temp > 40); length(ind)
[~,j] = find(temp2 > 40); length(j)
tempnum = fnum(j);
tmpnum = 0; k = 1;
for ii = 1:length(tempnum)
if tempnum(ii) ~= tmpnum
    tmpnum(k,1) = tempnum(ii); %#ok<SAGROW>
    k = k+1;
end
end
outofrangenum = [outofrangenum; tmpnum];

% for ii = 1:length(i)
%     figure(1); clf
%     plot(temp2(i(ii),:),-pres2(i(ii),:),'r.-'); hold on
%     plot(psal2(i(ii),:),-pres2(i(ii),:),'b.-');
%     temp2(i(ii),1:30)
%     pause
% end

ind = find(temp <= -3); length(ind)
[~,j] = find(temp2 <= -3); length(j)
tempnum = fnum(j);
tmpnum = 0; k = 1;
for ii = 1:length(tempnum)
if tempnum(ii) ~= tmpnum
    tmpnum(k,1) = tempnum(ii); %#ok<SAGROW>
    k = k+1;
end
end
outofrangenum = [outofrangenum; tmpnum];

% for ii = 1:length(i)
%     figure(1); clf
%     plot(temp2(i(ii),:),-pres2(i(ii),:),'r.-'); hold on
%     plot(psal2(i(ii),:),-pres2(i(ii),:),'b.-');
%     temp2(i(ii),1:30)
%     pause
% end

% look for out-of-range salinities
figure(1); clf; imagesc(psal); colorbar;
figure(2); clf; imagesc(psal2); colorbar

ind = find(psal > 39); length(ind)
[~,j] = find(psal2 > 39); length(j)
tempnum = fnum(j);
psalnum = 0; k = 1;
for ii = 1:length(tempnum)
if tempnum(ii) ~= psalnum
    psalnum(k,1) = tempnum(ii); %#ok<SAGROW>
    k = k+1;
end
end
outofrangenum = [outofrangenum; psalnum];

% for ii = 1:length(i)
%     figure(1); clf
% %     plot(temp2(i(ii),:),-pres2(i(ii),:),'r.-'); hold on
%     plot(psal2(:,j(ii)),-pres2(:,j(ii)),'b.-'); grid on
% %     psal2(i(ii),1:30)
%     pause
% end

ind = find(psal <= 20); length(ind)
[~,j] = find(psal2 <= 20); length(j)
tempnum = fnum(j);
psalnum = 0; k = 1;
for ii = 1:length(tempnum)
if tempnum(ii) ~= psalnum
    psalnum(k,1) = tempnum(ii); %#ok<SAGROW>
    k = k+1;
end
end
outofrangenum = [outofrangenum; psalnum];

% for ii = 1:length(i)
%     figure(1); clf
% %     plot(temp2(i(ii),:),-pres2(i(ii),:),'r.-'); hold on
%     plot(psal2(:,j(ii)),-pres2(:,j(ii)),'b.-'); grid on
% %     psal2(i(ii),1:30)
%     pause
% end

badnum = sort(outofrangenum);

clear psalnum tmpnum tempnum presnum i ii ind j k outofrangenum NaNnum ans

badnum = badnum(6:end)

fnumb = []; timeb = [];

badnum(i)
ind = find(fnum == badnum(i));
figure(1); clf; plot(temp2(:,ind),-pres2(:,ind),'.-'); grid on;
figure(2); clf; plot(psal2(:,ind),-pres2(:,ind),'.-'); grid on;
tmpt = temp2(:,ind); tmps = psal2(:,ind); tmpp = pres2(:,ind); tmpj = datestr(time(ind));
if sum(badnum(i) == grayfloat) > 0
  disp('greylisted');
  ii = find(badnum(i) == grayfloat);
  disp(param(ii)); disp(errtype(ii)); disp(start_date(ii)); disp(end_date(ii));
end

i = 1;
% set bad values to Nan
temp2(3,ind(8)) = NaN; 
psal2(isnan(temp2)) = NaN;
fnumb = [fnumb; fnum(ind(8))];
timeb = [timeb; time(ind(8))];

num_lowS = [fnumb timeb];
save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_lowS -APPEND

clear num_lowS fnumb timeb indb ind ans i ii tmpj tmpp tmps tmpt

% for ii = 1:length(ind)
%     ii
%     figure(1); clf; plot(psal2(ind,:)',-pres2(ind,:)','-','color',[.9 .9 .9]);
%     hold on; plot(psal2(ind(ii),:)',-pres2(ind(ii),:)','r.-'); grid on;
%     figure(2); clf;
%     axesm('mercator','MapLatLimit',[min(latitude(ind))-5 max(latitude(ind))+5],'MapLonLimit',[min(longitude(ind))-5 max(longitude(ind))+5])
%     plotm(latitude(ind),longitude(ind),'.-'); hold on;
%     plotm(latitude(ind(ii)),longitude(ind(ii)),'r.');
%     geoshow('landareas.shp','FaceColor',[.9 .9 .9])
%     setm(gca, 'grid','on','meridianlabel','on','parallellabel','on','plinelocation',10,'mlinelocation',10,'plabellocation',10,'mlabellocation',10)
%     axis tight;
%     figure(3); clf; plot(temp2(ind,:)',-pres2(ind,:)','-','color',[.9 .9 .9]);
%     hold on; plot(temp2(ind(ii),:)',-pres2(ind(ii),:)','r.-'); grid on;
%     pause
% end

% delete profiles that are all NaNs
clear ind indb; k = 1; l = 1;
for i = 1:length(fnum)
    t = temp2(:,i); s = psal2(:,i); p = pres2(:,i);
    ss = sum(~isnan(s)); tt = sum(~isnan(t)); pp = sum(~isnan(p));
    
    if ~(tt == 0 || ss == 0 || pp == 0)
        ind(k) = i; k = k+1; %#ok<SAGROW>
    else
        indb(l) = i; l = l+1; %#ok<SAGROW>
    end
end
clear t s p ss tt pp k i l
ind = ind'; indb = indb';

temp2 = temp2(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres2 = pres2(:,ind);
psal2 = psal2(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
pqc = pqc(:,ind); pres = pres(:,ind); psal = psal(:,ind); sqc = sqc(:,ind); temp = temp(:,ind); tqc = tqc(:,ind);
pres_adj = pres_adj(:,ind); pres_err = pres_err(:,ind); pres_qc = pres_qc(:,ind);
psal_adj = psal_adj(:,ind); psal_err = psal_err(:,ind); psal_qc = psal_qc(:,ind);
temp_adj = temp_adj(:,ind); temp_err = temp_err(:,ind); temp_qc = temp_qc(:,ind);
datamode = datamode(ind); time_adj = time_adj(ind);
clear ind indb 

%% find profiles with repeat pressures
clear ind1; k = 1;
for i = 1:length(fnum)
    p = pres2(:,i);
    for j = 1:length(p)
        if ~isnan(p(j))
            if sum(p == p(j)) > 1
                ind1(k) = i; %#ok<SAGROW>
                k = k+1;
            end
        end
    end
end
clear p i j k
ind1 = ind1';

% individual float numbers of profiles with repeat pressures
k = 1; ind = 0;
for i = 1:length(ind1)
    if fnum(ind1(i)) ~= ind
        ind(k) = fnum(ind1(i)); %#ok<SAGROW>
        k = k+1;
    end
end
ind = ind'; clear i k
length(ind)

for i = 1:length(ind)
    ii = find(fnum == ind(i));
    plot(temp2(:,ii),-pres2(:,ii),'.-');
    ind(i)
    pres2(1:2,ii)
    pause
end
clear ii i

i = 1;
temp2(71:72,ii(29)) = NaN; 
temp2(69,ii(32)) = NaN; 
psal2(isnan(temp2)) = NaN;

i = 2;
temp2(14,ii(10)) = NaN; 
psal2(isnan(temp2)) = NaN;

i = 3;
temp2(31,ii(18)) = NaN; 
psal2(isnan(temp2)) = NaN;

i = 4; 
i = 5; 
i = i+1

ind(i)
ii = find(fnum == ind(i));
figure(1); clf; plot(temp2(:,ii),-pres2(:,ii),'.-'); grid on;
figure(2); clf; plot(psal2(:,ii),-pres2(:,ii),'.-'); grid on;
tmpt = temp2(:,ii); tmps = psal2(:,ii); tmpp = pres2(:,ii); tmpj = datestr(time(ii));
if sum(ind(i) == grayfloat) > 0
  disp('greylisted');
  jj = find(ind(i) == grayfloat);
  disp(param(jj)); disp(errtype(jj)); disp(start_date(jj)); disp(end_date(jj));
end

clear ind2 ind3;
k = 1;
for iii = 1:length(ii)
p = pres2(:,ii(iii));
    for j = 1:length(p)
        if ~isnan(p(j))
            if sum(p == p(j)) > 1
                ind2(k) = iii;  %#ok<SAGROW>
                ind3(k) = j; %#ok<SAGROW>
                k = k+1;
            end
        end
    end
end
disp(ind2)

temp2(ind3(2),ii(ind2(1))) = NaN;
temp2(ind3(1:2),ii(ind2(1))) = NaN;
temp2(1:5,ii(ind2(1))) = NaN;
temp2(1:3,ii(11)) = NaN;
temp2(1:3,ii(15)) = NaN;
temp2(5,ii(16)) = NaN;
temp2(1:11,ii(18)) = NaN;

ind = find(isnan(pres2));
temp2(ind) = NaN;
psal2(ind) = NaN;

ind = find(isnan(psal2));
temp2(ind) = NaN;
pres2(ind) = NaN;

ind = find(isnan(temp2));
pres2(ind) = NaN;
psal2(ind) = NaN;

% find profiles with increasing pressures 
k = 1; ind1 = [];
for i = 1:length(fnum)
    p = pres2(:,i);
    p = p(~isnan(p));
    pd = p(2:end)-p(1:end-1);
    if sum(pd < 0) > 0
        ind1(k) = i; %#ok<SAGROW>
        k = k+1;
    end
end
clear ind i k p pd s t j
ind1 = ind1';

% individual float numbers of profiles with increasing pressures
k = 1; ind = 0;
for i = 1:length(ind1)
    if fnum(ind1(i)) ~= ind
        ind(k) = fnum(ind1(i)); %#ok<SAGROW>
        k = k+1;
    end
end
ind = ind'; clear i k
length(ind)

for i = 1:length(ind)
    ii = find(fnum == ind(i));
    plot(temp2(:,ii),-pres2(:,ii),'.-');
    ind(i)
    pres2(1:2,ii)
    pause
end
clear ii i

i = 1;
temp2(4,ii(ind2(1)))= NaN;

i = 2;
temp2(3,ii(ind2(1)))= NaN;

i = 3;
temp2(23,ii(ind2(1)))= NaN;

i = 4;
temp2(3,ii(ind2(1)))= NaN;

i = 5;
temp2(3,ii(ind2(1)))= NaN;

i = 6;
temp2(14,ii(ind2(1)))= NaN;

i = 7;
temp2(3,ii(ind2(1)))= NaN;

i = 8;
temp2(8,ii(ind2(1)))= NaN;

i = 9;
temp2(3,ii(ind2(1)))= NaN;

i = 10;
temp2(10,ii(ind2(1)))= NaN;

i = i+1
ind(i)
ii = find(fnum == ind(i));
figure(1); clf; plot(temp2(:,ii),-pres2(:,ii),'.-'); grid on;
figure(2); clf; plot(psal2(:,ii),-pres2(:,ii),'.-'); grid on;
tmpt = temp2(:,ii); tmps = psal2(:,ii); tmpp = pres2(:,ii); tmpj = datestr(time(ii));
if sum(ind(i) == grayfloat) > 0
  disp('greylisted');
  jj = find(ind(i) == grayfloat);
  disp(param(jj)); disp(errtype(jj)); disp(start_date(jj)); disp(end_date(jj));
end

k = 1; clear ind2
for iii = 1:length(ii)
    p = pres2(:,ii(iii));
    p = p(~isnan(p));
    pd = p(2:end)-p(1:end-1);
    if sum(pd < 0) > 0
        ind2(k) = iii; %#ok<SAGROW>
        k = k+1;
    end
end
disp(ind2)

ind = find(isnan(pres2));
temp2(ind) = NaN;
psal2(ind) = NaN;

ind = find(isnan(psal2));
temp2(ind) = NaN;
pres2(ind) = NaN;

ind = find(isnan(temp2));
pres2(ind) = NaN;
psal2(ind) = NaN;

clear ans i ii ind ind1 jj badnum tmpp tmpt tmpj tmps

%% plot all profiles and look for errant ones
figure(1); clf
plot(temp2,-pres2,'.-'); grid on; axis tight;

figure(2); clf 
plot(psal2,-pres2,'.-'); grid on; axis tight;

[i,j] = find(pres2 > 1900 & temp2 > 3);
pres2(i,j) = NaN;

[i,j] = find(pres2 > 1000 & psal2 > 34.75);
pres2(i,j) = NaN;

% i = 6;
% temp2(:,ind(12)) = NaN;
% psal2(isnan(temp2)) = NaN;
% fnumb = [fnumb; fnum(ind(12))]; timeb = [timeb; time(ind(12))];
% 
% badnum(i)
% ind = find(fnum == badnum(i));
% figure(1); clf; plot(temp2(:,ind),-pres2(:,ind),'.-'); grid on;
% figure(2); clf; plot(psal2(:,ind),-pres2(:,ind),'.-'); grid on;
% tmpt = temp2(:,ind); tmps = psal2(:,ind); tmpp = pres2(:,ind); tmpj = datestr(time(ind));
% if sum(badnum(i) == grayfloat) > 0
%   disp('greylisted');
%   ii = find(badnum(i) == grayfloat);
%   disp(param(ii)); disp(errtype(ii)); disp(start_date(ii)); disp(end_date(ii));
% end
% 
% num_errant = [fnumb timeb];
% save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_errant -APPEND
% clear num_errant fnumb timeb i ii ind indnum j jj badnum ans tmpj tmpp tmps tmpt


% need t, s, and p at every point - make NaN's consistent across all data
ind = find(isnan(pres2));
temp2(ind) = NaN;
psal2(ind) = NaN;

ind = find(isnan(psal2));
temp2(ind) = NaN;
pres2(ind) = NaN;

ind = find(isnan(temp2));
pres2(ind) = NaN;
psal2(ind) = NaN;

clear ind

% delete profiles that are all NaNs
clear ind indb; k = 1; l = 1;
for i = 1:length(fnum)
    t = temp2(:,i); s = psal2(:,i); p = pres2(:,i);
    ss = sum(~isnan(s)); tt = sum(~isnan(t)); pp = sum(~isnan(p));
    
    if ~(tt == 0 || ss == 0 || pp == 0)
        ind(k) = i; k = k+1; %#ok<SAGROW>
    else
        indb(l) = i; l = l+1; %#ok<SAGROW>
    end
end
clear t s p ss tt pp k i l
ind = ind'; indb = indb';

temp2 = temp2(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres2 = pres2(:,ind);
psal2 = psal2(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
pqc = pqc(:,ind); pres = pres(:,ind); psal = psal(:,ind); sqc = sqc(:,ind); temp = temp(:,ind); tqc = tqc(:,ind);
pres_adj = pres_adj(:,ind); pres_err = pres_err(:,ind); pres_qc = pres_qc(:,ind);
psal_adj = psal_adj(:,ind); psal_err = psal_err(:,ind); psal_qc = psal_qc(:,ind);
temp_adj = temp_adj(:,ind); temp_err = temp_err(:,ind); temp_qc = temp_qc(:,ind);
datamode = datamode(ind); time_adj = time_adj(ind);
clear ind indb

clear ans badnum end_date errtype grayfloat ii ind ind1 j jj kk num_* param
clear tmpj tmpp tmps tmpt start_date
% 
%% maximum pressure
maxp = (max(abs(pres2),[],1))';

% find & remove profiles with maxp < 800 db
indp = find(maxp < 800);
per800 = 100*length(indp)/length(maxp)

num_800_3 = [fnum(indp) time(indp)];

ind = find(maxp >= 800);
temp2 = temp2(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres2 = pres2(:,ind);
psal2 = psal2(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
pqc = pqc(:,ind); pres = pres(:,ind); psal = psal(:,ind); sqc = sqc(:,ind); temp = temp(:,ind); tqc = tqc(:,ind);
pres_adj = pres_adj(:,ind); pres_err = pres_err(:,ind); pres_qc = pres_qc(:,ind);
psal_adj = psal_adj(:,ind); psal_err = psal_err(:,ind); psal_qc = psal_qc(:,ind);
temp_adj = temp_adj(:,ind); temp_err = temp_err(:,ind); temp_qc = temp_qc(:,ind);
datamode = datamode(ind); time_adj = time_adj(ind);

clear ans ind indp per800 t2 

% % minimum pressure
% minp = (min(abs(pres2),[],1))';
% 
% % find & remove profiles with maxp < 800 db
% indp = find(minp > 400);
% per400 = 100*length(indp)/length(minp)
% 
% num_400 = [fnum(indp) time(indp)];
% 
% ind = find(minp <= 400);
% temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
% psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
% temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
% position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
% prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
% temp2 = temp2(:,ind); psal2 = psal2(:,ind); pres2 = pres2(:,ind);
% tqc = tqc(:,ind); pqc = pqc(:,ind); sqc = sqc(:,ind);
% 
% clear ans ind indp per800 t2 minp maxp per400

save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_800_3 -APPEND
clear num_800_3 num_400 maxp

clear ind3 jjj outofrangenum psalnum presnum tempnum tmpnum

%% remove fsi sensors
ind = find((type == 842 | (type == 847 | (type == 852 | type == 857))));

k = 1; indnum = 0;
for i = 1:length(fnum(ind))
    if fnum(ind(i)) ~= indnum
        indnum(k) = fnum(ind(i)); %#ok<SAGROW>
        k = k+1;
    end
end
indnum = indnum';
clear i k

cd /home/alison/Matlab/M' files'/Global' Mapping 2'/

for i = 1:length(indnum)
    fn = num2str(indnum(i));
    
    filename = strcat('/home/alison/Data/Argo Profiles/Float/North Pacific/','Profiles/',fn,'_prof.nc');
    nc = netcdf.open(char(filename),'NC_NOWRITE');
    pi = ncchar(nc,'PI_NAME');
    PI(i,:) = (pi(:,1))';
    netcdf.close(nc);
end

indnum2 = NaN(length(indnum),1);
for i = 1:length(indnum)
   
    if strcmp(PI(i,1:5),'BRECK')
        indnum2(i) = 1;
    else indnum2(i) = 0;
    end
end

indnum = indnum(indnum2==0);

for i = 1:length(indnum)

    ind = find(fnum == indnum(i));
    time(ind) = NaN;
end

ind = find(~isnan(time));

temp2 = temp2(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres2 = pres2(:,ind);
psal2 = psal2(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
pqc = pqc(:,ind); pres = pres(:,ind); psal = psal(:,ind); sqc = sqc(:,ind); temp = temp(:,ind); tqc = tqc(:,ind);
pres_adj = pres_adj(:,ind); pres_err = pres_err(:,ind); pres_qc = pres_qc(:,ind);
psal_adj = psal_adj(:,ind); psal_err = psal_err(:,ind); psal_qc = psal_qc(:,ind);
temp_adj = temp_adj(:,ind); temp_err = temp_err(:,ind); temp_qc = temp_qc(:,ind);
datamode = datamode(ind); time_adj = time_adj(ind);

clear PI ans filename fn i ind indnum indnum2 nc pi j


%% %%%%%%%%%%%%%%%
save('/home/alison/Matlab/Mat files/Global Mapping 2/North Pacific/prof_2005_qc_all.mat');   

clear juld_qc percentNaN position_qc pqc pres pres_qc prof_* psal psal_qc sqc dir_name indnum temp temp_qc tqc
pres = pres2; psal = psal2; temp = temp2; clear pres2 psal2 temp2 pres_adj pres_qc_adj psal_adj psal_qc_adj temp_adj temp_qc_adj time_adj

save('/home/alison/Matlab/Mat files/Global Mapping 2/North Pacific/prof_2005_qc.mat');

clear; load /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat
save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats_2005.mat



%% quality control

% 2006
clear; pack
load /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/prof_2006.mat;
load /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/prof_adj_2006.mat;
load /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/prof_err_2006.mat;

pres_qc = pres_qc';
psal_qc = psal_qc';
temp_qc = temp_qc';
pres_qc_adj = pres_qc_adj';
psal_qc_adj = psal_qc_adj';
temp_qc_adj = temp_qc_adj';

longitude(longitude < 0) = longitude(longitude < 0) + 360;

% plot locations of data
figure(1); clf
axesm('mercator','MapLatLimit',[-5 80],'MapLonLimit',[100 280])
plotm(latitude,longitude,'.');
geoshow('landareas.shp','FaceColor',[.9 .9 .9])
setm(gca, 'grid','on','meridianlabel','on','parallellabel','on','plinelocation',5,'mlinelocation',5,'plabellocation',5,'mlabellocation',5)
axis tight

figure(2); clf
plot(longitude,latitude,'.'); grid on;

% lat/lon range
ind = find(longitude <= 290 & longitude >= 100 & latitude <= 65 & latitude >= 0);
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

% remove gulf of mexico & atlantic profiles
ind = find(~(longitude >= 270 & latitude >= 16));
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

% remove sea of japan profiles
ind = find(~((longitude < 138 & latitude < 44.1 & latitude > 34.5) | (longitude < 141 & latitude < 44.1 & latitude > 37))); 
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

ind = find(~(longitude > 255 & latitude > 19 & latitude < 19.5));
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

ind = find(~(longitude < 126.5 & latitude > 9.05 & latitude < 9.4));
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

ind = find(~(longitude < 120.5 & latitude > 18 & latitude < 22.5));
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

ind = find(~(longitude < 125.625 & latitude > 11.315 & latitude < 11.33));
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);


% plot locations of data
figure(1); clf
axesm('mercator','MapLatLimit',[-5 80],'MapLonLimit',[100 280])
plotm(latitude,longitude,'.');
geoshow('landareas.shp','FaceColor',[.9 .9 .9])
setm(gca, 'grid','on','meridianlabel','on','parallellabel','on','plinelocation',10,'mlinelocation',10,'plabellocation',10,'mlabellocation',10)
axis tight

% look at maximum pressures
maxp = (max(abs(pres),[],1))';

% find & remove profiles with maxp < 800 db
indp = find(maxp < 800);
per800 = 100*length(indp)/length(maxp) %#ok<*NOPTS,*NASGU>

num_800 = [fnum(indp) time(indp)];

ind = find(maxp >= 800);
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

clear ans ind indp per800 t2 

save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_800
clear num_800 maxp

% delete all Nan profiles
clear ind; k = 1; 
for i = 1:length(fnum)
    t = temp(:,i); s = psal(:,i); p = pres(:,i);
    ss = sum(~isnan(s)); tt = sum(~isnan(t)); pp = sum(~isnan(p));
    
    if ~(tt == 0 || ss == 0 || pp == 0)
        ind(k) = i; k = k+1; %#ok<NAGROW>
    end
end
clear t s p ss tt pp k i
ind = ind';
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

indt = find(sum(~isnan(temp),2) == 0); %#ok<MXFND>
indp = find(sum(~isnan(pres),2) == 0); %#ok<MXFND>
inds = find(sum(~isnan(psal),2) == 0); %#ok<MXFND>
ind = max([min(indt) min(indp) min(inds)]);
temp = temp(1:ind,:); psal = psal(1:ind,:); pres = pres(1:ind,:);
temp_qc = temp_qc(1:ind,:); pres_qc = pres_qc(1:ind,:); psal_qc = psal_qc(1:ind,:);
temp_adj = temp_adj(1:ind,:); psal_adj = psal_adj(1:ind,:); pres_adj = pres_adj(1:ind,:);
temp_qc_adj = temp_qc_adj(1:ind,:); pres_qc_adj = pres_qc_adj(1:ind,:); psal_qc_adj = psal_qc_adj(1:ind,:);
temp_err = temp_err(1:ind,:); psal_err = psal_err(1:ind,:); pres_err = pres_err(1:ind,:);
clear indt indp inds ind

save /home/alison/Matlab/M'at files'/Global' Mapping 2'/North' Pacific'/tmp.mat

% graylist floats
cd('/home/alison/Matlab/M files/Global Mapping 2/');
[grayfloat,~,start_date,end_date,errtype] = graylist('/home/alison/Matlab/Mat files/Global Mapping 2/North Pacific/');


%% %%%% quality control flags %%%%

% examine flags on temperature data
[mm,nn] = size(temp);
tqc = NaN(size(temp));
for j = 1:nn
    for i = 1:mm
        if ~isnan(temp(i,j))
            if (strcmp(temp_qc_adj(i,j),'1') || strcmp(temp_qc_adj(i,j),'0') || strcmp(temp_qc_adj(i,j),'2') || strcmp(temp_qc_adj(i,j),' ') || double(temp_qc_adj(i,j))==0)
                tqc(i,j) = 1;
            else
                tqc(i,j) = 0;
            end
        end
    end
end

% examine flags on salinity data
[mm,nn] = size(psal);
sqc = NaN(size(psal));
for j = 1:nn
    for i = 1:mm
        if ~isnan(psal(i,j))
            if (strcmp(psal_qc_adj(i,j),'1') || strcmp(psal_qc_adj(i,j),'0') || strcmp(psal_qc_adj(i,j),'2') || strcmp(psal_qc_adj(i,j),' ') || double(psal_qc_adj(i,j))==0)
                sqc(i,j) = 1;
            else
                sqc(i,j) = 0;
            end
        end
    end
end

% examine flags on pressure data
[mm,nn] = size(pres);
pqc = NaN(size(pres));
for j = 1:nn
    for i = 1:mm
        if ~isnan(pres(i,j))
            if (strcmp(pres_qc_adj(i,j),'1') || strcmp(pres_qc_adj(i,j),'0') || strcmp(pres_qc_adj(i,j),'2') || strcmp(pres_qc_adj(i,j),' ') || double(pres_qc_adj(i,j))==0)
                pqc(i,j) = 1;
            else
                pqc(i,j) = 0;
            end
        end
    end
end

clear mm nn i j


% display number of profiles with one or more bad qc flag
num_bad_temp = length(find(sum(tqc == 0,1) > 0))
num_bad_pres = length(find(sum(pqc == 0,1) > 0))
num_bad_psal = length(find(sum(sqc == 0,1) > 0))
clear num_bad_*

% make 2nd variable set with data with qc ~= 0,1,2,' ' set to NaN
temp2 = temp_adj; temp2(tqc == 0) = NaN;
psal2 = psal_adj; psal2(sqc == 0) = NaN;
pres2 = pres_adj; pres2(pqc == 0) = NaN;

% need t, s, and p at every point - make NaN's consistent across all data
ind = find(isnan(pres2));
temp2(ind) = NaN;
psal2(ind) = NaN;

ind = find(isnan(psal2));
temp2(ind) = NaN;
pres2(ind) = NaN;

ind = find(isnan(temp2));
pres2(ind) = NaN;
psal2(ind) = NaN;

% delete profiles that are all NaNs
clear ind indb; k = 1; l = 1;
for i = 1:length(fnum)
    t = temp2(:,i); s = psal2(:,i); p = pres2(:,i);
    ss = sum(~isnan(s)); tt = sum(~isnan(t)); pp = sum(~isnan(p));
    
    if ~(tt == 0 || ss == 0 || pp == 0)
        ind(k) = i; k = k+1; %#ok<SAGROW>
    else
        indb(l) = i; l = l+1; %#ok<SAGROW>
    end
end
clear t s p ss tt pp k i l
ind = ind'; indb = indb';

num_badQC = [fnum(indb) time(indb)];
save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_badQC -APPEND
clear indb num_badQC

temp2 = temp2(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres2 = pres2(:,ind);
psal2 = psal2(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
pqc = pqc(:,ind); pres = pres(:,ind); psal = psal(:,ind); sqc = sqc(:,ind); temp = temp(:,ind); tqc = tqc(:,ind);
pres_adj = pres_adj(:,ind); pres_err = pres_err(:,ind); pres_qc = pres_qc(:,ind);
psal_adj = psal_adj(:,ind); psal_err = psal_err(:,ind); psal_qc = psal_qc(:,ind);
temp_adj = temp_adj(:,ind); temp_err = temp_err(:,ind); temp_qc = temp_qc(:,ind);
datamode = datamode(ind); time_adj = time_adj(ind);

clear ind

% find & remove profiles with maxp < 800 db
maxp = (max(abs(pres2),[],1))';
indp = find(maxp < 800);
per800 = 100*length(indp)/length(maxp) %#ok<*NOPTS,*NASGU>

num_800_2 = [fnum(indp) time(indp)];

ind = find(maxp >= 800);
temp2 = temp2(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres2 = pres2(:,ind);
psal2 = psal2(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
pqc = pqc(:,ind); pres = pres(:,ind); psal = psal(:,ind); sqc = sqc(:,ind); temp = temp(:,ind); tqc = tqc(:,ind);
pres_adj = pres_adj(:,ind); pres_err = pres_err(:,ind); pres_qc = pres_qc(:,ind);
psal_adj = psal_adj(:,ind); psal_err = psal_err(:,ind); psal_qc = psal_qc(:,ind);
temp_adj = temp_adj(:,ind); temp_err = temp_err(:,ind); temp_qc = temp_qc(:,ind);
datamode = datamode(ind); time_adj = time_adj(ind);
clear ind

clear ans ind indp per800 t2 

save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_800_2
clear num_800_2 maxp

% % compute percentage of good values in each profile
% [~,nn] = size(pres2);
% percentNaN = NaN(nn,3);
% for i = 1:nn
%     p = pqc(:,i); t = tqc(:,i); s = sqc(:,i);
%     np = sum(p == 0); nt = sum(t == 0); ns = sum(s == 0); 
%     npp = sum(~isnan(p)); ntt = sum(~isnan(t)); nss = sum(~isnan(s));
%     percentNaN(i,:) = 100*[np./npp nt/ntt ns/nss];
% end
% clear i p np lastp mm nn npp t nt ntt s ns nss
% 
% ind = find(percentNaN(:,1) > 25 | percentNaN(:,2) > 25 | percentNaN(:,3) > 25);
% for ii = 1:length(ind)
%     fnum(ind(ii))
%     figure(1); clf
%     plot(temp(:,ind(ii)),-pres(:,ind(ii)),'m.-'); hold on;
%     plot(temp2(:,ind(ii)),-pres2(:,ind(ii)),'r.-');
%     disp(temp_qc(:,ind(ii))');
%     disp(pres_qc(:,ind(ii))');
%     disp(psal_qc(:,ind(ii))');
%     figure(2); clf
%     plot(psal(:,ind(ii)),-pres(:,ind(ii)),'c.-'); hold on;
%     plot(psal2(:,ind(ii)),-pres2(:,ind(ii)),'b.-');
%     pause
% end
% clear ii ans
% 
% num_percentNaN = [fnum(ind) time(ind)];
% save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_percentNaN -APPEND
% clear num_percentNaN
% 
% temp2(:,ind) = NaN; pres2(:,ind) = NaN; psal2(:,ind) = NaN;

% find & delete short profiles
n = sum(~isnan(temp2),1)';
% for i = 1:10
%     disp(i); ind = find(n == i); length(ind)
%     figure(1); clf
%     plot(temp2(:,ind),-pres2(:,ind),'.-');
%     pause
% end

ind = find(n <= 10 & n > 0); length(ind)
num_shortprof = [fnum(ind) time(ind)];
save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_shortprof -APPEND
clear num_shortprof

temp2(:,ind) = NaN; pres2(:,ind) = NaN; psal2(:,ind) = NaN;
clear n i ind ans

% delete profiles that are all NaNs
clear ind indb; k = 1; l = 1;
for i = 1:length(fnum)
    t = temp2(:,i); s = psal2(:,i); p = pres2(:,i);
    ss = sum(~isnan(s)); tt = sum(~isnan(t)); pp = sum(~isnan(p));
    
    if ~(tt == 0 || ss == 0 || pp == 0)
        ind(k) = i; k = k+1; %#ok<SAGROW>
    else
        indb(l) = i; l = l+1; %#ok<SAGROW>
    end
end
clear t s p ss tt pp k i l
ind = ind'; indb = indb';

temp2 = temp2(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres2 = pres2(:,ind);
psal2 = psal2(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
pqc = pqc(:,ind); pres = pres(:,ind); psal = psal(:,ind); sqc = sqc(:,ind); temp = temp(:,ind); tqc = tqc(:,ind);
pres_adj = pres_adj(:,ind); pres_err = pres_err(:,ind); pres_qc = pres_qc(:,ind);
psal_adj = psal_adj(:,ind); psal_err = psal_err(:,ind); psal_qc = psal_qc(:,ind);
temp_adj = temp_adj(:,ind); temp_err = temp_err(:,ind); temp_qc = temp_qc(:,ind);
datamode = datamode(ind); time_adj = time_adj(ind);
clear ind indb 

save /home/alison/Matlab/M'at files'/Global' Mapping 2'/North' Pacific'/tmp.mat

%% remove greylist floats
k = 1;
badfloats = 0; type_err = cell(1,1); bad_start = 0; bad_end = 0;
for i = 1:length(grayfloat)
    if sum(grayfloat(i) == fnum) > 0
    if grayfloat(i) ~= badfloats
        badfloats(k,1) = grayfloat(i); %#ok<SAGROW>
        type_err(k,1) = errtype(i);
        bad_start(k,1) = start_date(i); %#ok<SAGROW>
        bad_end(k,1) = end_date(i); %#ok<SAGROW>
        k = k+1;
    end
    end
end
clear k i

% ind = find(~strcmp(type_err,'potential problem with bin assignment'));
% badfloats = badfloats(ind);
% type_err = type_err(ind);
% bad_start = bad_start(ind);
% bad_end = bad_end(ind);

badstart = cellstr(num2str(bad_start));
badend = cellstr(num2str(bad_end));

for i = 1:length(badstart)
    bs = char(badstart(i));
    bs = datenum([sscanf(bs(1:4),'%f'),sscanf(bs(5:6),'%f'),sscanf(bs(7:8),'%f')]);
    bad_st(i,1) = bs; %#ok<SAGROW>
    be = char(badend(i));
    if ~strcmp(be,'     NaN') && ~strcmp(be,'NaN')
        be = datenum([sscanf(be(1:4),'%f'),sscanf(be(5:6),'%f'),sscanf(be(7:8),'%f')]);
        bad_end(i) = be;
    else
        bad_end(i) = NaN;
    end
end

fnumb = []; timeb = [];
for j = 1:length(badfloats)
    if isnan(bad_end(j))
        ind = find(~((fnum == badfloats(j))&(time > bad_st(j))));
        indb = find((fnum == badfloats(j)) & (time > bad_st(j)));
    else
        ind = find(~(fnum == badfloats(j) & time > bad_st(j) & time < bad_end(j)));
        indb = find(fnum == badfloats(j) & time > bad_st(j) & time < bad_end(j));
    end
fnumb = [fnumb; fnum(indb)]; %#ok<AGROW>
timeb = [timeb; time(indb)]; %#ok<AGROW>

temp2 = temp2(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres2 = pres2(:,ind);
psal2 = psal2(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
pqc = pqc(:,ind); pres = pres(:,ind); psal = psal(:,ind); sqc = sqc(:,ind); temp = temp(:,ind); tqc = tqc(:,ind);
pres_adj = pres_adj(:,ind); pres_err = pres_err(:,ind); pres_qc = pres_qc(:,ind);
psal_adj = psal_adj(:,ind); psal_err = psal_err(:,ind); psal_qc = psal_qc(:,ind);
temp_adj = temp_adj(:,ind); temp_err = temp_err(:,ind); temp_qc = temp_qc(:,ind);
datamode = datamode(ind); time_adj = time_adj(ind);

end

num_gray = [fnumb timeb];
save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_gray -APPEND

clear num_gray fnumb timeb indb ind
clear bad* end_date j ind param start_date be bs i 
clear ans errtype grayfloat ii percentNaN type_err

% graylist floats
cd('/home/alison/Matlab/M files/Global Mapping 2/');
[grayfloat,param,start_date,end_date,errtype] = graylist('/home/alison/Matlab/Mat files/Global Mapping 2/North Pacific/');

save /home/alison/Matlab/M'at files'/Global' Mapping 2'/North' Pacific'/tmp.mat

%% look for out-of-range pressures
figure(1); clf; imagesc(pres); colorbar;
figure(2); clf; imagesc(pres2); colorbar

outofrangenum = [];

ind = find(pres > 2200); length(ind)
[i,j] = find(pres2 > 2200); length(i)
tempnum = fnum(j);
presnum = 0; k = 1;
for ii = 1:length(tempnum)
if tempnum(ii) ~= presnum
    presnum(k,1) = tempnum(ii); %#ok<SAGROW>
    k = k+1;
end
end
outofrangenum = [outofrangenum; presnum];

% for ii = 1:length(i)
%     figure(1); clf
%     plot(temp2(i(ii),:),-pres2(i(ii),:),'r.-'); hold on
%     plot(psal2(i(ii),:),-pres2(i(ii),:),'b.-');
%     pres2(i(ii),1:30)
%     pause
% end

ind = find(pres < 0); length(ind)
[~,j] = find(pres2 < 0); length(j)
tempnum = fnum(j);
presnum = 0; k = 1;
for ii = 1:length(tempnum)
if tempnum(ii) ~= presnum
    presnum(k,1) = tempnum(ii); %#ok<SAGROW>
    k = k+1;
end
end
outofrangenum = [outofrangenum; presnum];

% for ii = 1:length(i)
%     figure(1); clf
%     plot(temp(i(ii),:),-pres(i(ii),:),'r.-'); hold on
%     plot(psal(i(ii),:),-pres(i(ii),:),'b.-');
%     pres2(i(ii),1:30)
%     pause
% end

% look for out-of-range temperatures
figure(1); clf; imagesc(temp); colorbar;
figure(2); clf; imagesc(temp2); colorbar

ind = find(temp > 40); length(ind)
[~,j] = find(temp2 > 40); length(j)
tempnum = fnum(j);
tmpnum = 0; k = 1;
for ii = 1:length(tempnum)
if tempnum(ii) ~= tmpnum
    tmpnum(k,1) = tempnum(ii); %#ok<SAGROW>
    k = k+1;
end
end
outofrangenum = [outofrangenum; tmpnum];

% for ii = 1:length(i)
%     figure(1); clf
%     plot(temp2(i(ii),:),-pres2(i(ii),:),'r.-'); hold on
%     plot(psal2(i(ii),:),-pres2(i(ii),:),'b.-');
%     temp2(i(ii),1:30)
%     pause
% end

ind = find(temp <= -3); length(ind)
[~,j] = find(temp2 <= -3); length(j)
tempnum = fnum(j);
tmpnum = 0; k = 1;
for ii = 1:length(tempnum)
if tempnum(ii) ~= tmpnum
    tmpnum(k,1) = tempnum(ii); %#ok<SAGROW>
    k = k+1;
end
end
outofrangenum = [outofrangenum; tmpnum];

% for ii = 1:length(i)
%     figure(1); clf
%     plot(temp2(i(ii),:),-pres2(i(ii),:),'r.-'); hold on
%     plot(psal2(i(ii),:),-pres2(i(ii),:),'b.-');
%     temp2(i(ii),1:30)
%     pause
% end

% look for out-of-range salinities
figure(1); clf; imagesc(psal); colorbar;
figure(2); clf; imagesc(psal2); colorbar

ind = find(psal > 39); length(ind)
[~,j] = find(psal2 > 39); length(j)
tempnum = fnum(j);
psalnum = 0; k = 1;
for ii = 1:length(tempnum)
if tempnum(ii) ~= psalnum
    psalnum(k,1) = tempnum(ii); %#ok<SAGROW>
    k = k+1;
end
end
outofrangenum = [outofrangenum; psalnum];

% for ii = 1:length(i)
%     figure(1); clf
% %     plot(temp2(i(ii),:),-pres2(i(ii),:),'r.-'); hold on
%     plot(psal2(:,j(ii)),-pres2(:,j(ii)),'b.-'); grid on
% %     psal2(i(ii),1:30)
%     pause
% end

ind = find(psal <= 20); length(ind)
[~,j] = find(psal2 <= 20); length(j)
tempnum = fnum(j);
psalnum = 0; k = 1;
for ii = 1:length(tempnum)
if tempnum(ii) ~= psalnum
    psalnum(k,1) = tempnum(ii); %#ok<SAGROW>
    k = k+1;
end
end
outofrangenum = [outofrangenum; psalnum];

% for ii = 1:length(i)
%     figure(1); clf
% %     plot(temp2(i(ii),:),-pres2(i(ii),:),'r.-'); hold on
%     plot(psal2(:,j(ii)),-pres2(:,j(ii)),'b.-'); grid on
% %     psal2(i(ii),1:30)
%     pause
% end

badnum = sort(outofrangenum);

clear psalnum tmpnum tempnum presnum i ii ind j k outofrangenum NaNnum ans

badnum = badnum(6:end)

fnumb = []; timeb = [];

badnum(i)
ind = find(fnum == badnum(i));
figure(1); clf; plot(temp2(:,ind),-pres2(:,ind),'.-'); grid on;
figure(2); clf; plot(psal2(:,ind),-pres2(:,ind),'.-'); grid on;
tmpt = temp2(:,ind); tmps = psal2(:,ind); tmpp = pres2(:,ind); tmpj = datestr(time(ind));
if sum(badnum(i) == grayfloat) > 0
  disp('greylisted');
  ii = find(badnum(i) == grayfloat);
  disp(param(ii)); disp(errtype(ii)); disp(start_date(ii)); disp(end_date(ii));
end

i = 1;
% set bad values to Nan
temp2(1,ind(28)) = NaN; 
psal2(isnan(temp2)) = NaN;
fnumb = [fnumb; fnum(ind(28))];
timeb = [timeb; time(ind(28))];

num_lowS = [fnumb timeb];
save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_lowS -APPEND

clear num_lowS fnumb timeb indb ind ans i ii tmpj tmpp tmps tmpt

% for ii = 1:length(ind)
%     ii
%     figure(1); clf; plot(psal2(ind,:)',-pres2(ind,:)','-','color',[.9 .9 .9]);
%     hold on; plot(psal2(ind(ii),:)',-pres2(ind(ii),:)','r.-'); grid on;
%     figure(2); clf;
%     axesm('mercator','MapLatLimit',[min(latitude(ind))-5 max(latitude(ind))+5],'MapLonLimit',[min(longitude(ind))-5 max(longitude(ind))+5])
%     plotm(latitude(ind),longitude(ind),'.-'); hold on;
%     plotm(latitude(ind(ii)),longitude(ind(ii)),'r.');
%     geoshow('landareas.shp','FaceColor',[.9 .9 .9])
%     setm(gca, 'grid','on','meridianlabel','on','parallellabel','on','plinelocation',10,'mlinelocation',10,'plabellocation',10,'mlabellocation',10)
%     axis tight;
%     figure(3); clf; plot(temp2(ind,:)',-pres2(ind,:)','-','color',[.9 .9 .9]);
%     hold on; plot(temp2(ind(ii),:)',-pres2(ind(ii),:)','r.-'); grid on;
%     pause
% end

% delete profiles that are all NaNs
clear ind indb; k = 1; l = 1;
for i = 1:length(fnum)
    t = temp2(:,i); s = psal2(:,i); p = pres2(:,i);
    ss = sum(~isnan(s)); tt = sum(~isnan(t)); pp = sum(~isnan(p));
    
    if ~(tt == 0 || ss == 0 || pp == 0)
        ind(k) = i; k = k+1; %#ok<SAGROW>
    else
        indb(l) = i; l = l+1; %#ok<SAGROW>
    end
end
clear t s p ss tt pp k i l
ind = ind'; indb = indb';

temp2 = temp2(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres2 = pres2(:,ind);
psal2 = psal2(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
pqc = pqc(:,ind); pres = pres(:,ind); psal = psal(:,ind); sqc = sqc(:,ind); temp = temp(:,ind); tqc = tqc(:,ind);
pres_adj = pres_adj(:,ind); pres_err = pres_err(:,ind); pres_qc = pres_qc(:,ind);
psal_adj = psal_adj(:,ind); psal_err = psal_err(:,ind); psal_qc = psal_qc(:,ind);
temp_adj = temp_adj(:,ind); temp_err = temp_err(:,ind); temp_qc = temp_qc(:,ind);
datamode = datamode(ind); time_adj = time_adj(ind);
clear ind indb 

%% find profiles with repeat pressures
clear ind1; k = 1;
for i = 1:length(fnum)
    p = pres2(:,i);
    for j = 1:length(p)
        if ~isnan(p(j))
            if sum(p == p(j)) > 1
                ind1(k) = i; %#ok<SAGROW>
                k = k+1;
            end
        end
    end
end
clear p i j k
ind1 = ind1';

% individual float numbers of profiles with repeat pressures
k = 1; ind = 0;
for i = 1:length(ind1)
    if fnum(ind1(i)) ~= ind
        ind(k) = fnum(ind1(i)); %#ok<SAGROW>
        k = k+1;
    end
end
ind = ind'; clear i k
length(ind)

% for i = 1:length(ind)
%     ii = find(fnum == ind(i));
%     plot(temp2(:,ii),-pres2(:,ii),'.-');
%     ind(i)
%     pres2(1:2,ii)
%     pause
% end
% clear ii i

i = 1;

i = i+1;
ind(i)
ii = find(fnum == ind(i));
figure(1); clf; plot(temp2(:,ii),-pres2(:,ii),'.-'); grid on;
figure(2); clf; plot(psal2(:,ii),-pres2(:,ii),'.-'); grid on;
tmpt = temp2(:,ii); tmps = psal2(:,ii); tmpp = pres2(:,ii); tmpj = datestr(time(ii));
if sum(ind(i) == grayfloat) > 0
  disp('greylisted');
  jj = find(ind(i) == grayfloat);
  disp(param(jj)); disp(errtype(jj)); disp(start_date(jj)); disp(end_date(jj));
end

clear ind2 ind3;
k = 1;
for iii = 1:length(ii)
p = pres2(:,ii(iii));
    for j = 1:length(p)
        if ~isnan(p(j))
            if sum(p == p(j)) > 1
                ind2(k) = iii;  %#ok<SAGROW>
                ind3(k) = j; %#ok<SAGROW>
                k = k+1;
            end
        end
    end
end
disp(ind2)

temp2(ind3(2:3),ii(ind2(1))) = NaN;
temp2(ind3(1),ii(ind2(1))) = NaN;
temp2(ind3(4),ii(ind2(3))) = NaN;
temp2(ind3(6),ii(ind2(3))) = NaN;
temp2(46:48,ii(ind2(3))) = NaN;

ind = find(isnan(pres2));
temp2(ind) = NaN;
psal2(ind) = NaN;

ind = find(isnan(psal2));
temp2(ind) = NaN;
pres2(ind) = NaN;

ind = find(isnan(temp2));
pres2(ind) = NaN;
psal2(ind) = NaN;

% find profiles with increasing pressures 
k = 1; ind1 = [];
for i = 1:length(fnum)
    p = pres2(:,i);
    p = p(~isnan(p));
    pd = p(2:end)-p(1:end-1);
    if sum(pd < 0) > 0
        ind1(k) = i; %#ok<SAGROW>
        k = k+1;
    end
end
clear ind i k p pd s t j
ind1 = ind1';

% individual float numbers of profiles with increasing pressures
k = 1; ind = 0;
for i = 1:length(ind1)
    if fnum(ind1(i)) ~= ind
        ind(k) = fnum(ind1(i)); %#ok<SAGROW>
        k = k+1;
    end
end
ind = ind'; clear i k
length(ind)

for i = 1:length(ind)
    ii = find(fnum == ind(i));
    plot(temp2(:,ii),-pres2(:,ii),'.-');
    ind(i)
    pres2(1:2,ii)
    pause
end
clear ii i

i = 1;
pres2(13,ii(33)) = NaN; 

i = 2;
pres2(8,ii(47)) = NaN; 

i = 3;
pres2(3,ii(13)) = NaN; 

i = 4;
pres2(45,ii(46)) = NaN; 

i = 5;
pres2(8,ii(43)) = NaN; 
pres2(3,ii(49)) = NaN; 

i = 6;
pres2(13,ii(33)) = NaN; 
pres2(3,ii(62)) = NaN; 

i = 7;
for j =1:length(ind2)
    pres2(52,ii(ind2(j))) = NaN; 
end

ind(i)
ii = find(fnum == ind(i));
figure(1); clf; plot(temp2(:,ii),-pres2(:,ii),'.-'); grid on;
figure(2); clf; plot(psal2(:,ii),-pres2(:,ii),'.-'); grid on;
tmpt = temp2(:,ii); tmps = psal2(:,ii); tmpp = pres2(:,ii); tmpj = datestr(time(ii));
if sum(ind(i) == grayfloat) > 0
  disp('greylisted');
  jj = find(ind(i) == grayfloat);
  disp(param(jj)); disp(errtype(jj)); disp(start_date(jj)); disp(end_date(jj));
end

k = 1; clear ind2
for iii = 1:length(ii)
    p = pres2(:,ii(iii));
    p = p(~isnan(p));
    pd = p(2:end)-p(1:end-1);
    if sum(pd < 0) > 0
        ind2(k) = iii; %#ok<SAGROW>
        k = k+1;
    end
end
disp(ind2)

ind = find(isnan(pres2));
temp2(ind) = NaN;
psal2(ind) = NaN;

ind = find(isnan(psal2));
temp2(ind) = NaN;
pres2(ind) = NaN;

ind = find(isnan(temp2));
pres2(ind) = NaN;
psal2(ind) = NaN;

clear ans i ii ind ind1 jj badnum tmpp tmpt tmpj tmps

%% plot all profiles and look for errant ones
figure(1); clf
plot(temp2,-pres2,'.-'); grid on; axis tight;

figure(2); clf 
plot(psal2,-pres2,'.-'); grid on; axis tight;

[i,j] = find(pres2 > 1350 & temp2 < 0);
pres2(i,j) = NaN;

[i,j] = find(pres2 > 1850 & temp2 < 3.5 & temp2 > 2.75);
for k = 1:length(i)
    pres2(i(k),j(k)) = NaN;
end

[i,j] = find(pres2 > 900 & temp2 < 10.6 & temp2 > 10.2);
for k = 1:length(i)
    pres2(i(k),j(k)) = NaN;
end

[i,j] = find(pres2 > 1480 & psal2 < 34.3);
for k = 1:length(i)
    pres2(i(k),j(k)) = NaN;
end

% i = 6;
% temp2(:,ind(12)) = NaN;
% psal2(isnan(temp2)) = NaN;
% fnumb = [fnumb; fnum(ind(12))]; timeb = [timeb; time(ind(12))];
% 
% badnum(i)
% ind = find(fnum == badnum(i));
% figure(1); clf; plot(temp2(:,ind),-pres2(:,ind),'.-'); grid on;
% figure(2); clf; plot(psal2(:,ind),-pres2(:,ind),'.-'); grid on;
% tmpt = temp2(:,ind); tmps = psal2(:,ind); tmpp = pres2(:,ind); tmpj = datestr(time(ind));
% if sum(badnum(i) == grayfloat) > 0
%   disp('greylisted');
%   ii = find(badnum(i) == grayfloat);
%   disp(param(ii)); disp(errtype(ii)); disp(start_date(ii)); disp(end_date(ii));
% end
% 
% num_errant = [fnumb timeb];
% save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_errant -APPEND
% clear num_errant fnumb timeb i ii ind indnum j jj badnum ans tmpj tmpp tmps tmpt


% need t, s, and p at every point - make NaN's consistent across all data
ind = find(isnan(pres2));
temp2(ind) = NaN;
psal2(ind) = NaN;

ind = find(isnan(psal2));
temp2(ind) = NaN;
pres2(ind) = NaN;

ind = find(isnan(temp2));
pres2(ind) = NaN;
psal2(ind) = NaN;

clear ind

% delete profiles that are all NaNs
clear ind indb; k = 1; l = 1;
for i = 1:length(fnum)
    t = temp2(:,i); s = psal2(:,i); p = pres2(:,i);
    ss = sum(~isnan(s)); tt = sum(~isnan(t)); pp = sum(~isnan(p));
    
    if ~(tt == 0 || ss == 0 || pp == 0)
        ind(k) = i; k = k+1; %#ok<SAGROW>
    else
        indb(l) = i; l = l+1; %#ok<SAGROW>
    end
end
clear t s p ss tt pp k i l
ind = ind'; indb = indb';

temp2 = temp2(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres2 = pres2(:,ind);
psal2 = psal2(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
pqc = pqc(:,ind); pres = pres(:,ind); psal = psal(:,ind); sqc = sqc(:,ind); temp = temp(:,ind); tqc = tqc(:,ind);
pres_adj = pres_adj(:,ind); pres_err = pres_err(:,ind); pres_qc = pres_qc(:,ind);
psal_adj = psal_adj(:,ind); psal_err = psal_err(:,ind); psal_qc = psal_qc(:,ind);
temp_adj = temp_adj(:,ind); temp_err = temp_err(:,ind); temp_qc = temp_qc(:,ind);
datamode = datamode(ind); time_adj = time_adj(ind);
clear ind indb

clear ans badnum end_date errtype grayfloat ii ind ind1 j jj kk num_* param
clear tmpj tmpp tmps tmpt start_date
% 
%% maximum pressure
maxp = (max(abs(pres2),[],1))';

% find & remove profiles with maxp < 800 db
indp = find(maxp < 800);
per800 = 100*length(indp)/length(maxp)

num_800_3 = [fnum(indp) time(indp)];

ind = find(maxp >= 800);
temp2 = temp2(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres2 = pres2(:,ind);
psal2 = psal2(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
pqc = pqc(:,ind); pres = pres(:,ind); psal = psal(:,ind); sqc = sqc(:,ind); temp = temp(:,ind); tqc = tqc(:,ind);
pres_adj = pres_adj(:,ind); pres_err = pres_err(:,ind); pres_qc = pres_qc(:,ind);
psal_adj = psal_adj(:,ind); psal_err = psal_err(:,ind); psal_qc = psal_qc(:,ind);
temp_adj = temp_adj(:,ind); temp_err = temp_err(:,ind); temp_qc = temp_qc(:,ind);
datamode = datamode(ind); time_adj = time_adj(ind);

clear ans ind indp per800 t2 

% % minimum pressure
% minp = (min(abs(pres2),[],1))';
% 
% % find & remove profiles with maxp < 800 db
% indp = find(minp > 400);
% per400 = 100*length(indp)/length(minp)
% 
% num_400 = [fnum(indp) time(indp)];
% 
% ind = find(minp <= 400);
% temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
% psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
% temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
% position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
% prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
% temp2 = temp2(:,ind); psal2 = psal2(:,ind); pres2 = pres2(:,ind);
% tqc = tqc(:,ind); pqc = pqc(:,ind); sqc = sqc(:,ind);
% 
% clear ans ind indp per800 t2 minp maxp per400

save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_800_3 -APPEND
clear num_800_3 num_400 maxp

clear ind3 jjj outofrangenum psalnum presnum tempnum tmpnum

%% remove fsi sensors
ind = find((type == 842 | (type == 847 | (type == 852 | type == 857))));

k = 1; indnum = 0;
for i = 1:length(fnum(ind))
    if fnum(ind(i)) ~= indnum
        indnum(k) = fnum(ind(i)); %#ok<SAGROW>
        k = k+1;
    end
end
indnum = indnum';
clear i k

cd /home/alison/Matlab/M' files'/Global' Mapping 2'/

for i = 1:length(indnum)
    fn = num2str(indnum(i));
    
    filename = strcat('/home/alison/Data/Argo Profiles/Float/North Pacific/','Profiles/',fn,'_prof.nc');
    nc = netcdf.open(char(filename),'NC_NOWRITE');
    pi = ncchar(nc,'PI_NAME');
    PI(i,:) = (pi(:,1))';
    netcdf.close(nc);
end

indnum2 = NaN(length(indnum),1);
for i = 1:length(indnum)
   
    if strcmp(PI(i,1:5),'BRECK')
        indnum2(i) = 1;
    else indnum2(i) = 0;
    end
end

indnum = indnum(indnum2==0);

for i = 1:length(indnum)

    ind = find(fnum == indnum(i));
    time(ind) = NaN;
end

ind = find(~isnan(time));

temp2 = temp2(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres2 = pres2(:,ind);
psal2 = psal2(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
pqc = pqc(:,ind); pres = pres(:,ind); psal = psal(:,ind); sqc = sqc(:,ind); temp = temp(:,ind); tqc = tqc(:,ind);
pres_adj = pres_adj(:,ind); pres_err = pres_err(:,ind); pres_qc = pres_qc(:,ind);
psal_adj = psal_adj(:,ind); psal_err = psal_err(:,ind); psal_qc = psal_qc(:,ind);
temp_adj = temp_adj(:,ind); temp_err = temp_err(:,ind); temp_qc = temp_qc(:,ind);
datamode = datamode(ind); time_adj = time_adj(ind);

clear PI ans filename fn i ind indnum indnum2 nc pi j


%% %%%%%%%%%%%%%%%
save('/home/alison/Matlab/Mat files/Global Mapping 2/North Pacific/prof_2006_qc_all.mat');   

clear juld_qc percentNaN position_qc pqc pres pres_qc prof_* psal psal_qc sqc dir_name indnum temp temp_qc tqc
pres = pres2; psal = psal2; temp = temp2; clear pres2 psal2 temp2 pres_adj pres_qc_adj psal_adj psal_qc_adj temp_adj temp_qc_adj time_adj

save('/home/alison/Matlab/Mat files/Global Mapping 2/North Pacific/prof_2006_qc.mat');

clear; load /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat
save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats_2006.mat


%% quality control

% 2007
clear; pack
load /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/prof_2007.mat;
load /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/prof_adj_2007.mat;
load /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/prof_err_2007.mat;

pres_qc = pres_qc';
psal_qc = psal_qc';
temp_qc = temp_qc';
pres_qc_adj = pres_qc_adj';
psal_qc_adj = psal_qc_adj';
temp_qc_adj = temp_qc_adj';

longitude(longitude < 0) = longitude(longitude < 0) + 360;

% plot locations of data
figure(1); clf
axesm('mercator','MapLatLimit',[-5 80],'MapLonLimit',[100 280])
plotm(latitude,longitude,'.');
geoshow('landareas.shp','FaceColor',[.9 .9 .9])
setm(gca, 'grid','on','meridianlabel','on','parallellabel','on','plinelocation',5,'mlinelocation',5,'plabellocation',5,'mlabellocation',5)
axis tight

figure(2); clf
plot(longitude,latitude,'.'); grid on;

% lat/lon range
ind = find(longitude <= 282 & longitude >= 120.5 & latitude <= 63 & latitude >= 0);
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

% remove sea of japan profiles
ind = find(~((longitude < 138 & latitude < 44.1 & latitude > 34.5) | (longitude < 141.1 & latitude < 44.1 & latitude > 37))); 
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

ind = find(~(longitude > 236 & latitude > 47 & latitude < 48));
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

ind = find(~(longitude < 126 & latitude > 3 & latitude < 9.6));
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

ind = find(~(longitude < 126.2 & latitude > 9.06 & latitude < 9.08));
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

ind = find(~(longitude < 125.625 & latitude > 11.31 & latitude < 11.33));
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);


% plot locations of data
figure(1); clf
axesm('mercator','MapLatLimit',[-5 80],'MapLonLimit',[100 280])
plotm(latitude,longitude,'.');
geoshow('landareas.shp','FaceColor',[.9 .9 .9])
setm(gca, 'grid','on','meridianlabel','on','parallellabel','on','plinelocation',10,'mlinelocation',10,'plabellocation',10,'mlabellocation',10)
axis tight

% look at maximum pressures
maxp = (max(abs(pres),[],1))';

% find & remove profiles with maxp < 800 db
indp = find(maxp < 800);
per800 = 100*length(indp)/length(maxp) %#ok<*NOPTS,*NASGU>

num_800 = [fnum(indp) time(indp)];

ind = find(maxp >= 800);
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

clear ans ind indp per800 t2 

save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_800
clear num_800 maxp

% delete all Nan profiles
clear ind; k = 1; 
for i = 1:length(fnum)
    t = temp(:,i); s = psal(:,i); p = pres(:,i);
    ss = sum(~isnan(s)); tt = sum(~isnan(t)); pp = sum(~isnan(p));
    
    if ~(tt == 0 || ss == 0 || pp == 0)
        ind(k) = i; k = k+1; %#ok<NAGROW>
    end
end
clear t s p ss tt pp k i
ind = ind';
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

indt = find(sum(~isnan(temp),2) == 0); %#ok<MXFND>
indp = find(sum(~isnan(pres),2) == 0); %#ok<MXFND>
inds = find(sum(~isnan(psal),2) == 0); %#ok<MXFND>
ind = max([min(indt) min(indp) min(inds)]);
temp = temp(1:ind,:); psal = psal(1:ind,:); pres = pres(1:ind,:);
temp_qc = temp_qc(1:ind,:); pres_qc = pres_qc(1:ind,:); psal_qc = psal_qc(1:ind,:);
temp_adj = temp_adj(1:ind,:); psal_adj = psal_adj(1:ind,:); pres_adj = pres_adj(1:ind,:);
temp_qc_adj = temp_qc_adj(1:ind,:); pres_qc_adj = pres_qc_adj(1:ind,:); psal_qc_adj = psal_qc_adj(1:ind,:);
temp_err = temp_err(1:ind,:); psal_err = psal_err(1:ind,:); pres_err = pres_err(1:ind,:);
clear indt indp inds ind

save /home/alison/Matlab/M'at files'/Global' Mapping 2'/North' Pacific'/tmp.mat

% graylist floats
cd('/home/alison/Matlab/M files/Global Mapping 2/');
[grayfloat,~,start_date,end_date,errtype] = graylist('/home/alison/Matlab/Mat files/Global Mapping 2/North Pacific/');


%% %%%% quality control flags %%%%

% examine flags on temperature data
[mm,nn] = size(temp);
tqc = NaN(size(temp));
for j = 1:nn
    for i = 1:mm
        if ~isnan(temp(i,j))
            if (strcmp(temp_qc_adj(i,j),'1') || strcmp(temp_qc_adj(i,j),'0') || strcmp(temp_qc_adj(i,j),'2') || strcmp(temp_qc_adj(i,j),' ') || double(temp_qc_adj(i,j))==0)
                tqc(i,j) = 1;
            else
                tqc(i,j) = 0;
            end
        end
    end
end

% examine flags on salinity data
[mm,nn] = size(psal);
sqc = NaN(size(psal));
for j = 1:nn
    for i = 1:mm
        if ~isnan(psal(i,j))
            if (strcmp(psal_qc_adj(i,j),'1') || strcmp(psal_qc_adj(i,j),'0') || strcmp(psal_qc_adj(i,j),'2') || strcmp(psal_qc_adj(i,j),' ') || double(psal_qc_adj(i,j))==0)
                sqc(i,j) = 1;
            else
                sqc(i,j) = 0;
            end
        end
    end
end

% examine flags on pressure data
[mm,nn] = size(pres);
pqc = NaN(size(pres));
for j = 1:nn
    for i = 1:mm
        if ~isnan(pres(i,j))
            if (strcmp(pres_qc_adj(i,j),'1') || strcmp(pres_qc_adj(i,j),'0') || strcmp(pres_qc_adj(i,j),'2') || strcmp(pres_qc_adj(i,j),' ') || double(pres_qc_adj(i,j))==0)
                pqc(i,j) = 1;
            else
                pqc(i,j) = 0;
            end
        end
    end
end

clear mm nn i j


% display number of profiles with one or more bad qc flag
num_bad_temp = length(find(sum(tqc == 0,1) > 0))
num_bad_pres = length(find(sum(pqc == 0,1) > 0))
num_bad_psal = length(find(sum(sqc == 0,1) > 0))
clear num_bad_*

% make 2nd variable set with data with qc ~= 0,1,2,' ' set to NaN
temp2 = temp_adj; temp2(tqc == 0) = NaN;
psal2 = psal_adj; psal2(sqc == 0) = NaN;
pres2 = pres_adj; pres2(pqc == 0) = NaN;

% need t, s, and p at every point - make NaN's consistent across all data
ind = find(isnan(pres2));
temp2(ind) = NaN;
psal2(ind) = NaN;

ind = find(isnan(psal2));
temp2(ind) = NaN;
pres2(ind) = NaN;

ind = find(isnan(temp2));
pres2(ind) = NaN;
psal2(ind) = NaN;

% delete profiles that are all NaNs
clear ind indb; k = 1; l = 1;
for i = 1:length(fnum)
    t = temp2(:,i); s = psal2(:,i); p = pres2(:,i);
    ss = sum(~isnan(s)); tt = sum(~isnan(t)); pp = sum(~isnan(p));
    
    if ~(tt == 0 || ss == 0 || pp == 0)
        ind(k) = i; k = k+1; %#ok<SAGROW>
    else
        indb(l) = i; l = l+1; %#ok<SAGROW>
    end
end
clear t s p ss tt pp k i l
ind = ind'; indb = indb';

num_badQC = [fnum(indb) time(indb)];
save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_badQC -APPEND
clear indb num_badQC

temp2 = temp2(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres2 = pres2(:,ind);
psal2 = psal2(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
pqc = pqc(:,ind); pres = pres(:,ind); psal = psal(:,ind); sqc = sqc(:,ind); temp = temp(:,ind); tqc = tqc(:,ind);
pres_adj = pres_adj(:,ind); pres_err = pres_err(:,ind); pres_qc = pres_qc(:,ind);
psal_adj = psal_adj(:,ind); psal_err = psal_err(:,ind); psal_qc = psal_qc(:,ind);
temp_adj = temp_adj(:,ind); temp_err = temp_err(:,ind); temp_qc = temp_qc(:,ind);
datamode = datamode(ind); time_adj = time_adj(ind);

clear ind

% find & remove profiles with maxp < 800 db
maxp = (max(abs(pres2),[],1))';
indp = find(maxp < 800);
per800 = 100*length(indp)/length(maxp) %#ok<*NOPTS,*NASGU>

num_800_2 = [fnum(indp) time(indp)];

ind = find(maxp >= 800);
temp2 = temp2(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres2 = pres2(:,ind);
psal2 = psal2(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
pqc = pqc(:,ind); pres = pres(:,ind); psal = psal(:,ind); sqc = sqc(:,ind); temp = temp(:,ind); tqc = tqc(:,ind);
pres_adj = pres_adj(:,ind); pres_err = pres_err(:,ind); pres_qc = pres_qc(:,ind);
psal_adj = psal_adj(:,ind); psal_err = psal_err(:,ind); psal_qc = psal_qc(:,ind);
temp_adj = temp_adj(:,ind); temp_err = temp_err(:,ind); temp_qc = temp_qc(:,ind);
datamode = datamode(ind); time_adj = time_adj(ind);
clear ind

clear ans ind indp per800 t2 

save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_800_2
clear num_800_2 maxp

% % compute percentage of good values in each profile
% [~,nn] = size(pres2);
% percentNaN = NaN(nn,3);
% for i = 1:nn
%     p = pqc(:,i); t = tqc(:,i); s = sqc(:,i);
%     np = sum(p == 0); nt = sum(t == 0); ns = sum(s == 0); 
%     npp = sum(~isnan(p)); ntt = sum(~isnan(t)); nss = sum(~isnan(s));
%     percentNaN(i,:) = 100*[np./npp nt/ntt ns/nss];
% end
% clear i p np lastp mm nn npp t nt ntt s ns nss
% 
% ind = find(percentNaN(:,1) > 25 | percentNaN(:,2) > 25 | percentNaN(:,3) > 25);
% for ii = 1:length(ind)
%     fnum(ind(ii))
%     figure(1); clf
%     plot(temp(:,ind(ii)),-pres(:,ind(ii)),'m.-'); hold on;
%     plot(temp2(:,ind(ii)),-pres2(:,ind(ii)),'r.-');
%     disp(temp_qc(:,ind(ii))');
%     disp(pres_qc(:,ind(ii))');
%     disp(psal_qc(:,ind(ii))');
%     figure(2); clf
%     plot(psal(:,ind(ii)),-pres(:,ind(ii)),'c.-'); hold on;
%     plot(psal2(:,ind(ii)),-pres2(:,ind(ii)),'b.-');
%     pause
% end
% clear ii ans
% 
% num_percentNaN = [fnum(ind) time(ind)];
% save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_percentNaN -APPEND
% clear num_percentNaN
% 
% temp2(:,ind) = NaN; pres2(:,ind) = NaN; psal2(:,ind) = NaN;

% find & delete short profiles
n = sum(~isnan(temp2),1)';
% for i = 1:10
%     disp(i); ind = find(n == i); length(ind)
%     figure(1); clf
%     plot(temp2(:,ind),-pres2(:,ind),'.-');
%     pause
% end

ind = find(n <= 10 & n > 0); length(ind)
num_shortprof = [fnum(ind) time(ind)];
save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_shortprof -APPEND
clear num_shortprof

temp2(:,ind) = NaN; pres2(:,ind) = NaN; psal2(:,ind) = NaN;
clear n i ind ans

% delete profiles that are all NaNs
clear ind indb; k = 1; l = 1;
for i = 1:length(fnum)
    t = temp2(:,i); s = psal2(:,i); p = pres2(:,i);
    ss = sum(~isnan(s)); tt = sum(~isnan(t)); pp = sum(~isnan(p));
    
    if ~(tt == 0 || ss == 0 || pp == 0)
        ind(k) = i; k = k+1; %#ok<SAGROW>
    else
        indb(l) = i; l = l+1; %#ok<SAGROW>
    end
end
clear t s p ss tt pp k i l
ind = ind'; indb = indb';

temp2 = temp2(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres2 = pres2(:,ind);
psal2 = psal2(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
pqc = pqc(:,ind); pres = pres(:,ind); psal = psal(:,ind); sqc = sqc(:,ind); temp = temp(:,ind); tqc = tqc(:,ind);
pres_adj = pres_adj(:,ind); pres_err = pres_err(:,ind); pres_qc = pres_qc(:,ind);
psal_adj = psal_adj(:,ind); psal_err = psal_err(:,ind); psal_qc = psal_qc(:,ind);
temp_adj = temp_adj(:,ind); temp_err = temp_err(:,ind); temp_qc = temp_qc(:,ind);
datamode = datamode(ind); time_adj = time_adj(ind);
clear ind indb 

save /home/alison/Matlab/M'at files'/Global' Mapping 2'/North' Pacific'/tmp.mat

%% remove greylist floats
k = 1;
badfloats = 0; type_err = cell(1,1); bad_start = 0; bad_end = 0;
for i = 1:length(grayfloat)
    if sum(grayfloat(i) == fnum) > 0
    if grayfloat(i) ~= badfloats
        badfloats(k,1) = grayfloat(i); %#ok<SAGROW>
        type_err(k,1) = errtype(i);
        bad_start(k,1) = start_date(i); %#ok<SAGROW>
        bad_end(k,1) = end_date(i); %#ok<SAGROW>
        k = k+1;
    end
    end
end
clear k i

% ind = find(~strcmp(type_err,'potential problem with bin assignment'));
% badfloats = badfloats(ind);
% type_err = type_err(ind);
% bad_start = bad_start(ind);
% bad_end = bad_end(ind);

badstart = cellstr(num2str(bad_start));
badend = cellstr(num2str(bad_end));

for i = 1:length(badstart)
    bs = char(badstart(i));
    bs = datenum([sscanf(bs(1:4),'%f'),sscanf(bs(5:6),'%f'),sscanf(bs(7:8),'%f')]);
    bad_st(i,1) = bs; %#ok<SAGROW>
    be = char(badend(i));
    if ~strcmp(be,'     NaN') && ~strcmp(be,'NaN')
        be = datenum([sscanf(be(1:4),'%f'),sscanf(be(5:6),'%f'),sscanf(be(7:8),'%f')]);
        bad_end(i) = be;
    else
        bad_end(i) = NaN;
    end
end

fnumb = []; timeb = [];
for j = 1:length(badfloats)
    if isnan(bad_end(j))
        ind = find(~((fnum == badfloats(j))&(time > bad_st(j))));
        indb = find((fnum == badfloats(j)) & (time > bad_st(j)));
    else
        ind = find(~(fnum == badfloats(j) & time > bad_st(j) & time < bad_end(j)));
        indb = find(fnum == badfloats(j) & time > bad_st(j) & time < bad_end(j));
    end
fnumb = [fnumb; fnum(indb)]; %#ok<AGROW>
timeb = [timeb; time(indb)]; %#ok<AGROW>

temp2 = temp2(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres2 = pres2(:,ind);
psal2 = psal2(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
pqc = pqc(:,ind); pres = pres(:,ind); psal = psal(:,ind); sqc = sqc(:,ind); temp = temp(:,ind); tqc = tqc(:,ind);
pres_adj = pres_adj(:,ind); pres_err = pres_err(:,ind); pres_qc = pres_qc(:,ind);
psal_adj = psal_adj(:,ind); psal_err = psal_err(:,ind); psal_qc = psal_qc(:,ind);
temp_adj = temp_adj(:,ind); temp_err = temp_err(:,ind); temp_qc = temp_qc(:,ind);
datamode = datamode(ind); time_adj = time_adj(ind);

end

num_gray = [fnumb timeb];
save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_gray -APPEND

clear num_gray fnumb timeb indb ind
clear bad* end_date j ind param start_date be bs i 
clear ans errtype grayfloat ii percentNaN type_err

% graylist floats
cd('/home/alison/Matlab/M files/Global Mapping 2/');
[grayfloat,param,start_date,end_date,errtype] = graylist('/home/alison/Matlab/Mat files/Global Mapping 2/North Pacific/');

save /home/alison/Matlab/M'at files'/Global' Mapping 2'/North' Pacific'/tmp.mat

%% look for out-of-range pressures
figure(1); clf; imagesc(pres); colorbar;
figure(2); clf; imagesc(pres2); colorbar

outofrangenum = [];

ind = find(pres > 2200); length(ind)
[i,j] = find(pres2 > 2200); length(i)
tempnum = fnum(j);
presnum = 0; k = 1;
for ii = 1:length(tempnum)
if tempnum(ii) ~= presnum
    presnum(k,1) = tempnum(ii); %#ok<SAGROW>
    k = k+1;
end
end
outofrangenum = [outofrangenum; presnum];

% for ii = 1:length(i)
%     figure(1); clf
%     plot(temp2(i(ii),:),-pres2(i(ii),:),'r.-'); hold on
%     plot(psal2(i(ii),:),-pres2(i(ii),:),'b.-');
%     pres2(i(ii),1:30)
%     pause
% end

ind = find(pres < 0); length(ind)
[~,j] = find(pres2 < 0); length(j)
tempnum = fnum(j);
presnum = 0; k = 1;
for ii = 1:length(tempnum)
if tempnum(ii) ~= presnum
    presnum(k,1) = tempnum(ii); %#ok<SAGROW>
    k = k+1;
end
end
outofrangenum = [outofrangenum; presnum];

% for ii = 1:length(i)
%     figure(1); clf
%     plot(temp(i(ii),:),-pres(i(ii),:),'r.-'); hold on
%     plot(psal(i(ii),:),-pres(i(ii),:),'b.-');
%     pres2(i(ii),1:30)
%     pause
% end

% look for out-of-range temperatures
figure(1); clf; imagesc(temp); colorbar;
figure(2); clf; imagesc(temp2); colorbar

ind = find(temp > 40); length(ind)
[~,j] = find(temp2 > 40); length(j)
tempnum = fnum(j);
tmpnum = 0; k = 1;
for ii = 1:length(tempnum)
if tempnum(ii) ~= tmpnum
    tmpnum(k,1) = tempnum(ii); %#ok<SAGROW>
    k = k+1;
end
end
outofrangenum = [outofrangenum; tmpnum];

% for ii = 1:length(i)
%     figure(1); clf
%     plot(temp2(i(ii),:),-pres2(i(ii),:),'r.-'); hold on
%     plot(psal2(i(ii),:),-pres2(i(ii),:),'b.-');
%     temp2(i(ii),1:30)
%     pause
% end

ind = find(temp <= -3); length(ind)
[~,j] = find(temp2 <= -3); length(j)
tempnum = fnum(j);
tmpnum = 0; k = 1;
for ii = 1:length(tempnum)
if tempnum(ii) ~= tmpnum
    tmpnum(k,1) = tempnum(ii); %#ok<SAGROW>
    k = k+1;
end
end
outofrangenum = [outofrangenum; tmpnum];

% for ii = 1:length(i)
%     figure(1); clf
%     plot(temp2(i(ii),:),-pres2(i(ii),:),'r.-'); hold on
%     plot(psal2(i(ii),:),-pres2(i(ii),:),'b.-');
%     temp2(i(ii),1:30)
%     pause
% end

% look for out-of-range salinities
figure(1); clf; imagesc(psal); colorbar;
figure(2); clf; imagesc(psal2); colorbar

ind = find(psal > 39); length(ind)
[~,j] = find(psal2 > 39); length(j)
tempnum = fnum(j);
psalnum = 0; k = 1;
for ii = 1:length(tempnum)
if tempnum(ii) ~= psalnum
    psalnum(k,1) = tempnum(ii); %#ok<SAGROW>
    k = k+1;
end
end
outofrangenum = [outofrangenum; psalnum];

% for ii = 1:length(i)
%     figure(1); clf
% %     plot(temp2(i(ii),:),-pres2(i(ii),:),'r.-'); hold on
%     plot(psal2(:,j(ii)),-pres2(:,j(ii)),'b.-'); grid on
% %     psal2(i(ii),1:30)
%     pause
% end

ind = find(psal <= 20); length(ind)
[~,j] = find(psal2 <= 20); length(j)
tempnum = fnum(j);
psalnum = 0; k = 1;
for ii = 1:length(tempnum)
if tempnum(ii) ~= psalnum
    psalnum(k,1) = tempnum(ii); %#ok<SAGROW>
    k = k+1;
end
end
outofrangenum = [outofrangenum; psalnum];

% for ii = 1:length(i)
%     figure(1); clf
% %     plot(temp2(i(ii),:),-pres2(i(ii),:),'r.-'); hold on
%     plot(psal2(:,j(ii)),-pres2(:,j(ii)),'b.-'); grid on
% %     psal2(i(ii),1:30)
%     pause
% end

badnum = sort(outofrangenum);

clear psalnum tmpnum tempnum presnum i ii ind j k outofrangenum NaNnum ans

badnum = badnum(6:end)

fnumb = []; timeb = [];

badnum(i)
ind = find(fnum == badnum(i));
figure(1); clf; plot(temp2(:,ind),-pres2(:,ind),'.-'); grid on;
figure(2); clf; plot(psal2(:,ind),-pres2(:,ind),'.-'); grid on;
tmpt = temp2(:,ind); tmps = psal2(:,ind); tmpp = pres2(:,ind); tmpj = datestr(time(ind));
if sum(badnum(i) == grayfloat) > 0
  disp('greylisted');
  ii = find(badnum(i) == grayfloat);
  disp(param(ii)); disp(errtype(ii)); disp(start_date(ii)); disp(end_date(ii));
end

i = 1;
% set bad values to Nan
temp2(1,ind(24)) = NaN; 
psal2(isnan(temp2)) = NaN;
fnumb = [fnumb; fnum(ind(24))];
timeb = [timeb; time(ind(24))];

i = 2;
% set bad values to Nan
temp2(1,ind(26)) = NaN; 
temp2(2:3,ind(35)) = NaN;
psal2(isnan(temp2)) = NaN;
fnumb = [fnumb; fnum(ind(26))];
timeb = [timeb; time(ind(26))];

num_lowS = [fnumb timeb];
save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_lowS -APPEND

clear num_lowS fnumb timeb indb ind ans i ii tmpj tmpp tmps tmpt

% for ii = 1:length(ind)
%     ii
%     figure(1); clf; plot(psal2(ind,:)',-pres2(ind,:)','-','color',[.9 .9 .9]);
%     hold on; plot(psal2(ind(ii),:)',-pres2(ind(ii),:)','r.-'); grid on;
%     figure(2); clf;
%     axesm('mercator','MapLatLimit',[min(latitude(ind))-5 max(latitude(ind))+5],'MapLonLimit',[min(longitude(ind))-5 max(longitude(ind))+5])
%     plotm(latitude(ind),longitude(ind),'.-'); hold on;
%     plotm(latitude(ind(ii)),longitude(ind(ii)),'r.');
%     geoshow('landareas.shp','FaceColor',[.9 .9 .9])
%     setm(gca, 'grid','on','meridianlabel','on','parallellabel','on','plinelocation',10,'mlinelocation',10,'plabellocation',10,'mlabellocation',10)
%     axis tight;
%     figure(3); clf; plot(temp2(ind,:)',-pres2(ind,:)','-','color',[.9 .9 .9]);
%     hold on; plot(temp2(ind(ii),:)',-pres2(ind(ii),:)','r.-'); grid on;
%     pause
% end

% delete profiles that are all NaNs
clear ind indb; k = 1; l = 1;
for i = 1:length(fnum)
    t = temp2(:,i); s = psal2(:,i); p = pres2(:,i);
    ss = sum(~isnan(s)); tt = sum(~isnan(t)); pp = sum(~isnan(p));
    
    if ~(tt == 0 || ss == 0 || pp == 0)
        ind(k) = i; k = k+1; %#ok<SAGROW>
    else
        indb(l) = i; l = l+1; %#ok<SAGROW>
    end
end
clear t s p ss tt pp k i l
ind = ind'; indb = indb';

temp2 = temp2(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres2 = pres2(:,ind);
psal2 = psal2(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
pqc = pqc(:,ind); pres = pres(:,ind); psal = psal(:,ind); sqc = sqc(:,ind); temp = temp(:,ind); tqc = tqc(:,ind);
pres_adj = pres_adj(:,ind); pres_err = pres_err(:,ind); pres_qc = pres_qc(:,ind);
psal_adj = psal_adj(:,ind); psal_err = psal_err(:,ind); psal_qc = psal_qc(:,ind);
temp_adj = temp_adj(:,ind); temp_err = temp_err(:,ind); temp_qc = temp_qc(:,ind);
datamode = datamode(ind); time_adj = time_adj(ind);
clear ind indb 

%% find profiles with repeat pressures
clear ind1; k = 1;
for i = 1:length(fnum)
    p = pres2(:,i);
    for j = 1:length(p)
        if ~isnan(p(j))
            if sum(p == p(j)) > 1
                ind1(k) = i; %#ok<SAGROW>
                k = k+1;
            end
        end
    end
end
clear p i j k
ind1 = ind1';

% individual float numbers of profiles with repeat pressures
k = 1; ind = 0;
for i = 1:length(ind1)
    if fnum(ind1(i)) ~= ind
        ind(k) = fnum(ind1(i)); %#ok<SAGROW>
        k = k+1;
    end
end
ind = ind'; clear i k
length(ind)

% for i = 1:length(ind)
%     ii = find(fnum == ind(i));
%     plot(temp2(:,ii),-pres2(:,ii),'.-');
%     ind(i)
%     pres2(1:2,ii)
%     pause
% end
% clear ii i

i = 1;

i = i+1;
ind(i)
ii = find(fnum == ind(i));
figure(1); clf; plot(temp2(:,ii),-pres2(:,ii),'.-'); grid on;
figure(2); clf; plot(psal2(:,ii),-pres2(:,ii),'.-'); grid on;
tmpt = temp2(:,ii); tmps = psal2(:,ii); tmpp = pres2(:,ii); tmpj = datestr(time(ii));
if sum(ind(i) == grayfloat) > 0
  disp('greylisted');
  jj = find(ind(i) == grayfloat);
  disp(param(jj)); disp(errtype(jj)); disp(start_date(jj)); disp(end_date(jj));
end

clear ind2 ind3;
k = 1;
for iii = 1:length(ii)
p = pres2(:,ii(iii));
    for j = 1:length(p)
        if ~isnan(p(j))
            if sum(p == p(j)) > 1
                ind2(k) = iii;  %#ok<SAGROW>
                ind3(k) = j; %#ok<SAGROW>
                k = k+1;
            end
        end
    end
end
disp(ind2)

temp2(ind3(2:3),ii(ind2(1))) = NaN;
temp2(ind3(1),ii(ind2(1))) = NaN;
temp2(ind3(4),ii(ind2(3))) = NaN;
temp2(ind3(6),ii(ind2(3))) = NaN;
temp2(46:48,ii(ind2(3))) = NaN;

ind = find(isnan(pres2));
temp2(ind) = NaN;
psal2(ind) = NaN;

ind = find(isnan(psal2));
temp2(ind) = NaN;
pres2(ind) = NaN;

ind = find(isnan(temp2));
pres2(ind) = NaN;
psal2(ind) = NaN;

% find profiles with increasing pressures 
k = 1; ind1 = [];
for i = 1:length(fnum)
    p = pres2(:,i);
    p = p(~isnan(p));
    pd = p(2:end)-p(1:end-1);
    if sum(pd < 0) > 0
        ind1(k) = i; %#ok<SAGROW>
        k = k+1;
    end
end
clear ind i k p pd s t j
ind1 = ind1';

% individual float numbers of profiles with increasing pressures
k = 1; ind = 0;
for i = 1:length(ind1)
    if fnum(ind1(i)) ~= ind
        ind(k) = fnum(ind1(i)); %#ok<SAGROW>
        k = k+1;
    end
end
ind = ind'; clear i k
length(ind)

for i = 1:length(ind)
    ii = find(fnum == ind(i));
    plot(temp2(:,ii),-pres2(:,ii),'.-');
    ind(i)
    pres2(1:2,ii)
    pause
end
clear ii i

i = 1;
pres2(63,ii(18)) = NaN; 

i = 2;
pres2(41,ii(54)) = NaN; 
pres2(9,ii(54)) = NaN; 

i = 3;
pres2(44,ii(41)) = NaN; 

i = 4;
pres2(25,ii(2)) = NaN; 
pres2(3,ii(9)) = NaN; 
pres2(1,ii(20)) = NaN; 

i = 5;
pres2(33,ii(2)) = NaN; 

i = 6;
for j =1:length(ind2)
    pres2(52,ii(ind2(j))) = NaN; 
end

i = 7;
pres2(6,ii(21)) = NaN; 


ind(i)
ii = find(fnum == ind(i));
figure(1); clf; plot(temp2(:,ii),-pres2(:,ii),'.-'); grid on;
figure(2); clf; plot(psal2(:,ii),-pres2(:,ii),'.-'); grid on;
tmpt = temp2(:,ii); tmps = psal2(:,ii); tmpp = pres2(:,ii); tmpj = datestr(time(ii));
if sum(ind(i) == grayfloat) > 0
  disp('greylisted');
  jj = find(ind(i) == grayfloat);
  disp(param(jj)); disp(errtype(jj)); disp(start_date(jj)); disp(end_date(jj));
end

k = 1; clear ind2
for iii = 1:length(ii)
    p = pres2(:,ii(iii));
    p = p(~isnan(p));
    pd = p(2:end)-p(1:end-1);
    if sum(pd < 0) > 0
        ind2(k) = iii; %#ok<SAGROW>
        k = k+1;
    end
end
disp(ind2)

ind = find(isnan(pres2));
temp2(ind) = NaN;
psal2(ind) = NaN;

ind = find(isnan(psal2));
temp2(ind) = NaN;
pres2(ind) = NaN;

ind = find(isnan(temp2));
pres2(ind) = NaN;
psal2(ind) = NaN;

clear ans i ii ind ind1 jj badnum tmpp tmpt tmpj tmps
clear iii ind2

%% plot all profiles and look for errant ones
figure(1); clf
plot(temp2,-pres2,'.-'); grid on; axis tight;

figure(2); clf 
plot(psal2,-pres2,'.-'); grid on; axis tight;

[i,j] = find(pres2 > 1920 & temp2 > 2.75);
for k = 1:length(i)
    pres2(i(k),j(k)) = NaN;
end

[i,j] = find(pres2 > 1500 & psal2 < 34.4);
for k = 1:length(i)
    pres2(i(k),j(k)) = NaN;
end


% i = 6;
% temp2(:,ind(12)) = NaN;
% psal2(isnan(temp2)) = NaN;
% fnumb = [fnumb; fnum(ind(12))]; timeb = [timeb; time(ind(12))];
% 
% badnum(i)
% ind = find(fnum == badnum(i));
% figure(1); clf; plot(temp2(:,ind),-pres2(:,ind),'.-'); grid on;
% figure(2); clf; plot(psal2(:,ind),-pres2(:,ind),'.-'); grid on;
% tmpt = temp2(:,ind); tmps = psal2(:,ind); tmpp = pres2(:,ind); tmpj = datestr(time(ind));
% if sum(badnum(i) == grayfloat) > 0
%   disp('greylisted');
%   ii = find(badnum(i) == grayfloat);
%   disp(param(ii)); disp(errtype(ii)); disp(start_date(ii)); disp(end_date(ii));
% end
% 
% num_errant = [fnumb timeb];
% save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_errant -APPEND
% clear num_errant fnumb timeb i ii ind indnum j jj badnum ans tmpj tmpp tmps tmpt


% need t, s, and p at every point - make NaN's consistent across all data
ind = find(isnan(pres2));
temp2(ind) = NaN;
psal2(ind) = NaN;

ind = find(isnan(psal2));
temp2(ind) = NaN;
pres2(ind) = NaN;

ind = find(isnan(temp2));
pres2(ind) = NaN;
psal2(ind) = NaN;

clear ind

% delete profiles that are all NaNs
clear ind indb; k = 1; l = 1;
for i = 1:length(fnum)
    t = temp2(:,i); s = psal2(:,i); p = pres2(:,i);
    ss = sum(~isnan(s)); tt = sum(~isnan(t)); pp = sum(~isnan(p));
    
    if ~(tt == 0 || ss == 0 || pp == 0)
        ind(k) = i; k = k+1; %#ok<SAGROW>
    else
        indb(l) = i; l = l+1; %#ok<SAGROW>
    end
end
clear t s p ss tt pp k i l
ind = ind'; indb = indb';

temp2 = temp2(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres2 = pres2(:,ind);
psal2 = psal2(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
pqc = pqc(:,ind); pres = pres(:,ind); psal = psal(:,ind); sqc = sqc(:,ind); temp = temp(:,ind); tqc = tqc(:,ind);
pres_adj = pres_adj(:,ind); pres_err = pres_err(:,ind); pres_qc = pres_qc(:,ind);
psal_adj = psal_adj(:,ind); psal_err = psal_err(:,ind); psal_qc = psal_qc(:,ind);
temp_adj = temp_adj(:,ind); temp_err = temp_err(:,ind); temp_qc = temp_qc(:,ind);
datamode = datamode(ind); time_adj = time_adj(ind);
clear ind indb

clear ans badnum end_date errtype grayfloat ii ind ind1 j jj kk num_* param
clear tmpj tmpp tmps tmpt start_date
% 
%% maximum pressure
maxp = (max(abs(pres2),[],1))';

% find & remove profiles with maxp < 800 db
indp = find(maxp < 800);
per800 = 100*length(indp)/length(maxp)

num_800_3 = [fnum(indp) time(indp)];

ind = find(maxp >= 800);
temp2 = temp2(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres2 = pres2(:,ind);
psal2 = psal2(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
pqc = pqc(:,ind); pres = pres(:,ind); psal = psal(:,ind); sqc = sqc(:,ind); temp = temp(:,ind); tqc = tqc(:,ind);
pres_adj = pres_adj(:,ind); pres_err = pres_err(:,ind); pres_qc = pres_qc(:,ind);
psal_adj = psal_adj(:,ind); psal_err = psal_err(:,ind); psal_qc = psal_qc(:,ind);
temp_adj = temp_adj(:,ind); temp_err = temp_err(:,ind); temp_qc = temp_qc(:,ind);
datamode = datamode(ind); time_adj = time_adj(ind);

clear ans ind indp per800 t2 

% % minimum pressure
% minp = (min(abs(pres2),[],1))';
% 
% % find & remove profiles with maxp < 800 db
% indp = find(minp > 400);
% per400 = 100*length(indp)/length(minp)
% 
% num_400 = [fnum(indp) time(indp)];
% 
% ind = find(minp <= 400);
% temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
% psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
% temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
% position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
% prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
% temp2 = temp2(:,ind); psal2 = psal2(:,ind); pres2 = pres2(:,ind);
% tqc = tqc(:,ind); pqc = pqc(:,ind); sqc = sqc(:,ind);
% 
% clear ans ind indp per800 t2 minp maxp per400

save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_800_3 -APPEND
clear num_800_3 num_400 maxp

clear ind3 jjj outofrangenum psalnum presnum tempnum tmpnum

%% remove fsi sensors
ind = find((type == 842 | (type == 847 | (type == 852 | type == 857))));

k = 1; indnum = 0;
for i = 1:length(fnum(ind))
    if fnum(ind(i)) ~= indnum
        indnum(k) = fnum(ind(i)); %#ok<SAGROW>
        k = k+1;
    end
end
indnum = indnum';
clear i k

cd /home/alison/Matlab/M' files'/Global' Mapping 2'/

for i = 1:length(indnum)
    fn = num2str(indnum(i));
    
    filename = strcat('/home/alison/Data/Argo Profiles/Float/North Pacific/','Profiles/',fn,'_prof.nc');
    nc = netcdf.open(char(filename),'NC_NOWRITE');
    pi = ncchar(nc,'PI_NAME');
    PI(i,:) = (pi(:,1))';
    netcdf.close(nc);
end

indnum2 = NaN(length(indnum),1);
for i = 1:length(indnum)
   
    if strcmp(PI(i,1:5),'BRECK')
        indnum2(i) = 1;
    else indnum2(i) = 0;
    end
end

indnum = indnum(indnum2==0);

for i = 1:length(indnum)

    ind = find(fnum == indnum(i));
    time(ind) = NaN;
end

ind = find(~isnan(time));

temp2 = temp2(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres2 = pres2(:,ind);
psal2 = psal2(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
pqc = pqc(:,ind); pres = pres(:,ind); psal = psal(:,ind); sqc = sqc(:,ind); temp = temp(:,ind); tqc = tqc(:,ind);
pres_adj = pres_adj(:,ind); pres_err = pres_err(:,ind); pres_qc = pres_qc(:,ind);
psal_adj = psal_adj(:,ind); psal_err = psal_err(:,ind); psal_qc = psal_qc(:,ind);
temp_adj = temp_adj(:,ind); temp_err = temp_err(:,ind); temp_qc = temp_qc(:,ind);
datamode = datamode(ind); time_adj = time_adj(ind);

clear PI ans filename fn i ind indnum indnum2 nc pi j


%% %%%%%%%%%%%%%%%
save('/home/alison/Matlab/Mat files/Global Mapping 2/North Pacific/prof_2007_qc_all.mat');   

clear juld_qc percentNaN position_qc pqc pres pres_qc prof_* psal psal_qc sqc dir_name indnum temp temp_qc tqc
pres = pres2; psal = psal2; temp = temp2; clear pres2 psal2 temp2 pres_adj pres_qc_adj psal_adj psal_qc_adj temp_adj temp_qc_adj time_adj

save('/home/alison/Matlab/Mat files/Global Mapping 2/North Pacific/prof_2007_qc.mat');

clear; load /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat
save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats_2007.mat


%% quality control

% 2008
clear; pack
load /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/prof_2008.mat;
load /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/prof_adj_2008.mat;
load /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/prof_err_2008.mat;

pres_qc = pres_qc';
psal_qc = psal_qc';
temp_qc = temp_qc';
pres_qc_adj = pres_qc_adj';
psal_qc_adj = psal_qc_adj';
temp_qc_adj = temp_qc_adj';

longitude(longitude < 0) = longitude(longitude < 0) + 360;

% plot locations of data
figure(1); clf
axesm('mercator','MapLatLimit',[-5 80],'MapLonLimit',[100 280])
plotm(latitude,longitude,'.');
geoshow('landareas.shp','FaceColor',[.9 .9 .9])
setm(gca, 'grid','on','meridianlabel','on','parallellabel','on','plinelocation',5,'mlinelocation',5,'plabellocation',5,'mlabellocation',5)
axis tight

figure(2); clf
plot(longitude,latitude,'.'); grid on;

% lat/lon range
ind = find(longitude <= 282 & longitude >= 120.5 & latitude <= 63 & latitude >= 0);
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

% remove sea of japan profiles
ind = find(~((longitude < 138 & latitude < 44.1 & latitude > 34.5) | (longitude < 141.1 & latitude < 45.1 & latitude > 37))); 
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

ind = find(~(longitude > 236 & latitude > 47 & latitude < 48));
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

ind = find(~(longitude < 126.5 & latitude > 3 & latitude < 8));
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

ind = find(~(longitude < 123.5 & latitude > 9.06 & latitude < 9.5));
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

ind = find(~(longitude < 125.8 & latitude > 11.25 & latitude < 11.35));
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

ind = find(~(longitude < 126.2 & latitude > 9.06 & latitude < 9.08));
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

ind = find(~(longitude < 126.04 & latitude > 10.02 & latitude < 10.04));
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

ind = find(~(longitude > 243.49 & latitude > 31.4 & latitude < 31.42));
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

% plot locations of data
figure(1); clf
axesm('mercator','MapLatLimit',[-5 80],'MapLonLimit',[100 280])
plotm(latitude,longitude,'.');
geoshow('landareas.shp','FaceColor',[.9 .9 .9])
setm(gca, 'grid','on','meridianlabel','on','parallellabel','on','plinelocation',10,'mlinelocation',10,'plabellocation',10,'mlabellocation',10)
axis tight

% look at maximum pressures
maxp = (max(abs(pres),[],1))';

% find & remove profiles with maxp < 800 db
indp = find(maxp < 800);
per800 = 100*length(indp)/length(maxp) %#ok<*NOPTS,*NASGU>

num_800 = [fnum(indp) time(indp)];

ind = find(maxp >= 800);
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

clear ans ind indp per800 t2 

save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_800
clear num_800 maxp

% delete all Nan profiles
clear ind; k = 1; 
for i = 1:length(fnum)
    t = temp(:,i); s = psal(:,i); p = pres(:,i);
    ss = sum(~isnan(s)); tt = sum(~isnan(t)); pp = sum(~isnan(p));
    
    if ~(tt == 0 || ss == 0 || pp == 0)
        ind(k) = i; k = k+1; %#ok<NAGROW>
    end
end
clear t s p ss tt pp k i
ind = ind';
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

indt = find(sum(~isnan(temp),2) == 0); %#ok<MXFND>
indp = find(sum(~isnan(pres),2) == 0); %#ok<MXFND>
inds = find(sum(~isnan(psal),2) == 0); %#ok<MXFND>
ind = max([min(indt) min(indp) min(inds)]);
temp = temp(1:ind,:); psal = psal(1:ind,:); pres = pres(1:ind,:);
temp_qc = temp_qc(1:ind,:); pres_qc = pres_qc(1:ind,:); psal_qc = psal_qc(1:ind,:);
temp_adj = temp_adj(1:ind,:); psal_adj = psal_adj(1:ind,:); pres_adj = pres_adj(1:ind,:);
temp_qc_adj = temp_qc_adj(1:ind,:); pres_qc_adj = pres_qc_adj(1:ind,:); psal_qc_adj = psal_qc_adj(1:ind,:);
temp_err = temp_err(1:ind,:); psal_err = psal_err(1:ind,:); pres_err = pres_err(1:ind,:);
clear indt indp inds ind

save /home/alison/Matlab/M'at files'/Global' Mapping 2'/North' Pacific'/tmp.mat

% graylist floats
cd('/home/alison/Matlab/M files/Global Mapping 2/');
[grayfloat,~,start_date,end_date,errtype] = graylist('/home/alison/Matlab/Mat files/Global Mapping 2/North Pacific/');


%% %%%% quality control flags %%%%

% examine flags on temperature data
[mm,nn] = size(temp);
tqc = NaN(size(temp));
for j = 1:nn
    for i = 1:mm
        if ~isnan(temp(i,j))
            if (strcmp(temp_qc_adj(i,j),'1') || strcmp(temp_qc_adj(i,j),'0') || strcmp(temp_qc_adj(i,j),'2') || strcmp(temp_qc_adj(i,j),' ') || double(temp_qc_adj(i,j))==0)
                tqc(i,j) = 1;
            else
                tqc(i,j) = 0;
            end
        end
    end
end

% examine flags on salinity data
[mm,nn] = size(psal);
sqc = NaN(size(psal));
for j = 1:nn
    for i = 1:mm
        if ~isnan(psal(i,j))
            if (strcmp(psal_qc_adj(i,j),'1') || strcmp(psal_qc_adj(i,j),'0') || strcmp(psal_qc_adj(i,j),'2') || strcmp(psal_qc_adj(i,j),' ') || double(psal_qc_adj(i,j))==0)
                sqc(i,j) = 1;
            else
                sqc(i,j) = 0;
            end
        end
    end
end

% examine flags on pressure data
[mm,nn] = size(pres);
pqc = NaN(size(pres));
for j = 1:nn
    for i = 1:mm
        if ~isnan(pres(i,j))
            if (strcmp(pres_qc_adj(i,j),'1') || strcmp(pres_qc_adj(i,j),'0') || strcmp(pres_qc_adj(i,j),'2') || strcmp(pres_qc_adj(i,j),' ') || double(pres_qc_adj(i,j))==0)
                pqc(i,j) = 1;
            else
                pqc(i,j) = 0;
            end
        end
    end
end

clear mm nn i j


% display number of profiles with one or more bad qc flag
num_bad_temp = length(find(sum(tqc == 0,1) > 0))
num_bad_pres = length(find(sum(pqc == 0,1) > 0))
num_bad_psal = length(find(sum(sqc == 0,1) > 0))
clear num_bad_*

% make 2nd variable set with data with qc ~= 0,1,2,' ' set to NaN
temp2 = temp_adj; temp2(tqc == 0) = NaN;
psal2 = psal_adj; psal2(sqc == 0) = NaN;
pres2 = pres_adj; pres2(pqc == 0) = NaN;

% need t, s, and p at every point - make NaN's consistent across all data
ind = find(isnan(pres2));
temp2(ind) = NaN;
psal2(ind) = NaN;

ind = find(isnan(psal2));
temp2(ind) = NaN;
pres2(ind) = NaN;

ind = find(isnan(temp2));
pres2(ind) = NaN;
psal2(ind) = NaN;

% delete profiles that are all NaNs
clear ind indb; k = 1; l = 1;
for i = 1:length(fnum)
    t = temp2(:,i); s = psal2(:,i); p = pres2(:,i);
    ss = sum(~isnan(s)); tt = sum(~isnan(t)); pp = sum(~isnan(p));
    
    if ~(tt == 0 || ss == 0 || pp == 0)
        ind(k) = i; k = k+1; %#ok<SAGROW>
    else
        indb(l) = i; l = l+1; %#ok<SAGROW>
    end
end
clear t s p ss tt pp k i l
ind = ind'; indb = indb';

num_badQC = [fnum(indb) time(indb)];
save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_badQC -APPEND
clear indb num_badQC

temp2 = temp2(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres2 = pres2(:,ind);
psal2 = psal2(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
pqc = pqc(:,ind); pres = pres(:,ind); psal = psal(:,ind); sqc = sqc(:,ind); temp = temp(:,ind); tqc = tqc(:,ind);
pres_adj = pres_adj(:,ind); pres_err = pres_err(:,ind); pres_qc = pres_qc(:,ind);
psal_adj = psal_adj(:,ind); psal_err = psal_err(:,ind); psal_qc = psal_qc(:,ind);
temp_adj = temp_adj(:,ind); temp_err = temp_err(:,ind); temp_qc = temp_qc(:,ind);
datamode = datamode(ind); time_adj = time_adj(ind);

clear ind

% find & remove profiles with maxp < 800 db
maxp = (max(abs(pres2),[],1))';
indp = find(maxp < 800);
per800 = 100*length(indp)/length(maxp) %#ok<*NOPTS,*NASGU>

num_800_2 = [fnum(indp) time(indp)];

ind = find(maxp >= 800);
temp2 = temp2(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres2 = pres2(:,ind);
psal2 = psal2(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
pqc = pqc(:,ind); pres = pres(:,ind); psal = psal(:,ind); sqc = sqc(:,ind); temp = temp(:,ind); tqc = tqc(:,ind);
pres_adj = pres_adj(:,ind); pres_err = pres_err(:,ind); pres_qc = pres_qc(:,ind);
psal_adj = psal_adj(:,ind); psal_err = psal_err(:,ind); psal_qc = psal_qc(:,ind);
temp_adj = temp_adj(:,ind); temp_err = temp_err(:,ind); temp_qc = temp_qc(:,ind);
datamode = datamode(ind); time_adj = time_adj(ind);
clear ind

clear ans ind indp per800 t2 

save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_800_2
clear num_800_2 maxp

% % compute percentage of good values in each profile
% [~,nn] = size(pres2);
% percentNaN = NaN(nn,3);
% for i = 1:nn
%     p = pqc(:,i); t = tqc(:,i); s = sqc(:,i);
%     np = sum(p == 0); nt = sum(t == 0); ns = sum(s == 0); 
%     npp = sum(~isnan(p)); ntt = sum(~isnan(t)); nss = sum(~isnan(s));
%     percentNaN(i,:) = 100*[np./npp nt/ntt ns/nss];
% end
% clear i p np lastp mm nn npp t nt ntt s ns nss
% 
% ind = find(percentNaN(:,1) > 25 | percentNaN(:,2) > 25 | percentNaN(:,3) > 25);
% for ii = 1:length(ind)
%     fnum(ind(ii))
%     figure(1); clf
%     plot(temp(:,ind(ii)),-pres(:,ind(ii)),'m.-'); hold on;
%     plot(temp2(:,ind(ii)),-pres2(:,ind(ii)),'r.-');
%     disp(temp_qc(:,ind(ii))');
%     disp(pres_qc(:,ind(ii))');
%     disp(psal_qc(:,ind(ii))');
%     figure(2); clf
%     plot(psal(:,ind(ii)),-pres(:,ind(ii)),'c.-'); hold on;
%     plot(psal2(:,ind(ii)),-pres2(:,ind(ii)),'b.-');
%     pause
% end
% clear ii ans
% 
% num_percentNaN = [fnum(ind) time(ind)];
% save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_percentNaN -APPEND
% clear num_percentNaN
% 
% temp2(:,ind) = NaN; pres2(:,ind) = NaN; psal2(:,ind) = NaN;

% find & delete short profiles
n = sum(~isnan(temp2),1)';
% for i = 1:10
%     disp(i); ind = find(n == i); length(ind)
%     figure(1); clf
%     plot(temp2(:,ind),-pres2(:,ind),'.-');
%     pause
% end

ind = find(n <= 10 & n > 0); length(ind)
num_shortprof = [fnum(ind) time(ind)];
save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_shortprof -APPEND
clear num_shortprof

temp2(:,ind) = NaN; pres2(:,ind) = NaN; psal2(:,ind) = NaN;
clear n i ind ans

% delete profiles that are all NaNs
clear ind indb; k = 1; l = 1;
for i = 1:length(fnum)
    t = temp2(:,i); s = psal2(:,i); p = pres2(:,i);
    ss = sum(~isnan(s)); tt = sum(~isnan(t)); pp = sum(~isnan(p));
    
    if ~(tt == 0 || ss == 0 || pp == 0)
        ind(k) = i; k = k+1; %#ok<SAGROW>
    else
        indb(l) = i; l = l+1; %#ok<SAGROW>
    end
end
clear t s p ss tt pp k i l
ind = ind'; indb = indb';

temp2 = temp2(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres2 = pres2(:,ind);
psal2 = psal2(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
pqc = pqc(:,ind); pres = pres(:,ind); psal = psal(:,ind); sqc = sqc(:,ind); temp = temp(:,ind); tqc = tqc(:,ind);
pres_adj = pres_adj(:,ind); pres_err = pres_err(:,ind); pres_qc = pres_qc(:,ind);
psal_adj = psal_adj(:,ind); psal_err = psal_err(:,ind); psal_qc = psal_qc(:,ind);
temp_adj = temp_adj(:,ind); temp_err = temp_err(:,ind); temp_qc = temp_qc(:,ind);
datamode = datamode(ind); time_adj = time_adj(ind);
clear ind indb 

save /home/alison/Matlab/M'at files'/Global' Mapping 2'/North' Pacific'/tmp.mat

%% remove greylist floats
k = 1;
badfloats = 0; type_err = cell(1,1); bad_start = 0; bad_end = 0;
for i = 1:length(grayfloat)
    if sum(grayfloat(i) == fnum) > 0
    if grayfloat(i) ~= badfloats
        badfloats(k,1) = grayfloat(i); %#ok<SAGROW>
        type_err(k,1) = errtype(i);
        bad_start(k,1) = start_date(i); %#ok<SAGROW>
        bad_end(k,1) = end_date(i); %#ok<SAGROW>
        k = k+1;
    end
    end
end
clear k i

% ind = find(~strcmp(type_err,'potential problem with bin assignment'));
% badfloats = badfloats(ind);
% type_err = type_err(ind);
% bad_start = bad_start(ind);
% bad_end = bad_end(ind);

badstart = cellstr(num2str(bad_start));
badend = cellstr(num2str(bad_end));

for i = 1:length(badstart)
    bs = char(badstart(i));
    bs = datenum([sscanf(bs(1:4),'%f'),sscanf(bs(5:6),'%f'),sscanf(bs(7:8),'%f')]);
    bad_st(i,1) = bs; %#ok<SAGROW>
    be = char(badend(i));
    if ~strcmp(be,'     NaN') && ~strcmp(be,'NaN')
        be = datenum([sscanf(be(1:4),'%f'),sscanf(be(5:6),'%f'),sscanf(be(7:8),'%f')]);
        bad_end(i) = be;
    else
        bad_end(i) = NaN;
    end
end

fnumb = []; timeb = [];
for j = 1:length(badfloats)
    if isnan(bad_end(j))
        ind = find(~((fnum == badfloats(j))&(time > bad_st(j))));
        indb = find((fnum == badfloats(j)) & (time > bad_st(j)));
    else
        ind = find(~(fnum == badfloats(j) & time > bad_st(j) & time < bad_end(j)));
        indb = find(fnum == badfloats(j) & time > bad_st(j) & time < bad_end(j));
    end
fnumb = [fnumb; fnum(indb)]; %#ok<AGROW>
timeb = [timeb; time(indb)]; %#ok<AGROW>

temp2 = temp2(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres2 = pres2(:,ind);
psal2 = psal2(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
pqc = pqc(:,ind); pres = pres(:,ind); psal = psal(:,ind); sqc = sqc(:,ind); temp = temp(:,ind); tqc = tqc(:,ind);
pres_adj = pres_adj(:,ind); pres_err = pres_err(:,ind); pres_qc = pres_qc(:,ind);
psal_adj = psal_adj(:,ind); psal_err = psal_err(:,ind); psal_qc = psal_qc(:,ind);
temp_adj = temp_adj(:,ind); temp_err = temp_err(:,ind); temp_qc = temp_qc(:,ind);
datamode = datamode(ind); time_adj = time_adj(ind);

end

num_gray = [fnumb timeb];
save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_gray -APPEND

clear num_gray fnumb timeb indb ind
clear bad* end_date j ind param start_date be bs i 
clear ans errtype grayfloat ii percentNaN type_err

% graylist floats
cd('/home/alison/Matlab/M files/Global Mapping 2/');
[grayfloat,param,start_date,end_date,errtype] = graylist('/home/alison/Matlab/Mat files/Global Mapping 2/North Pacific/');

save /home/alison/Matlab/M'at files'/Global' Mapping 2'/North' Pacific'/tmp.mat

%% look for out-of-range pressures
figure(1); clf; imagesc(pres); colorbar;
figure(2); clf; imagesc(pres2); colorbar

outofrangenum = [];

ind = find(pres > 2200); length(ind)
[i,j] = find(pres2 > 2200); length(i)
tempnum = fnum(j);
presnum = 0; k = 1;
for ii = 1:length(tempnum)
if tempnum(ii) ~= presnum
    presnum(k,1) = tempnum(ii); %#ok<SAGROW>
    k = k+1;
end
end
outofrangenum = [outofrangenum; presnum];

% for ii = 1:length(i)
%     figure(1); clf
%     plot(temp2(i(ii),:),-pres2(i(ii),:),'r.-'); hold on
%     plot(psal2(i(ii),:),-pres2(i(ii),:),'b.-');
%     pres2(i(ii),1:30)
%     pause
% end

ind = find(pres < 0); length(ind)
[~,j] = find(pres2 < 0); length(j)
tempnum = fnum(j);
presnum = 0; k = 1;
for ii = 1:length(tempnum)
if tempnum(ii) ~= presnum
    presnum(k,1) = tempnum(ii); %#ok<SAGROW>
    k = k+1;
end
end
outofrangenum = [outofrangenum; presnum];

% for ii = 1:length(i)
%     figure(1); clf
%     plot(temp(i(ii),:),-pres(i(ii),:),'r.-'); hold on
%     plot(psal(i(ii),:),-pres(i(ii),:),'b.-');
%     pres2(i(ii),1:30)
%     pause
% end

% look for out-of-range temperatures
figure(1); clf; imagesc(temp); colorbar;
figure(2); clf; imagesc(temp2); colorbar

ind = find(temp > 40); length(ind)
[~,j] = find(temp2 > 40); length(j)
tempnum = fnum(j);
tmpnum = 0; k = 1;
for ii = 1:length(tempnum)
if tempnum(ii) ~= tmpnum
    tmpnum(k,1) = tempnum(ii); %#ok<SAGROW>
    k = k+1;
end
end
outofrangenum = [outofrangenum; tmpnum];

% for ii = 1:length(i)
%     figure(1); clf
%     plot(temp2(i(ii),:),-pres2(i(ii),:),'r.-'); hold on
%     plot(psal2(i(ii),:),-pres2(i(ii),:),'b.-');
%     temp2(i(ii),1:30)
%     pause
% end

ind = find(temp <= -3); length(ind)
[~,j] = find(temp2 <= -3); length(j)
tempnum = fnum(j);
tmpnum = 0; k = 1;
for ii = 1:length(tempnum)
if tempnum(ii) ~= tmpnum
    tmpnum(k,1) = tempnum(ii); %#ok<SAGROW>
    k = k+1;
end
end
outofrangenum = [outofrangenum; tmpnum];

% for ii = 1:length(i)
%     figure(1); clf
%     plot(temp2(i(ii),:),-pres2(i(ii),:),'r.-'); hold on
%     plot(psal2(i(ii),:),-pres2(i(ii),:),'b.-');
%     temp2(i(ii),1:30)
%     pause
% end

% look for out-of-range salinities
figure(1); clf; imagesc(psal); colorbar;
figure(2); clf; imagesc(psal2); colorbar

ind = find(psal > 39); length(ind)
[~,j] = find(psal2 > 39); length(j)
tempnum = fnum(j);
psalnum = 0; k = 1;
for ii = 1:length(tempnum)
if tempnum(ii) ~= psalnum
    psalnum(k,1) = tempnum(ii); %#ok<SAGROW>
    k = k+1;
end
end
outofrangenum = [outofrangenum; psalnum];

% for ii = 1:length(i)
%     figure(1); clf
% %     plot(temp2(i(ii),:),-pres2(i(ii),:),'r.-'); hold on
%     plot(psal2(:,j(ii)),-pres2(:,j(ii)),'b.-'); grid on
% %     psal2(i(ii),1:30)
%     pause
% end

ind = find(psal <= 20); length(ind)
[~,j] = find(psal2 <= 20); length(j)
tempnum = fnum(j);
psalnum = 0; k = 1;
for ii = 1:length(tempnum)
if tempnum(ii) ~= psalnum
    psalnum(k,1) = tempnum(ii); %#ok<SAGROW>
    k = k+1;
end
end
outofrangenum = [outofrangenum; psalnum];

% for ii = 1:length(i)
%     figure(1); clf
% %     plot(temp2(i(ii),:),-pres2(i(ii),:),'r.-'); hold on
%     plot(psal2(:,j(ii)),-pres2(:,j(ii)),'b.-'); grid on
% %     psal2(i(ii),1:30)
%     pause
% end

badnum = sort(outofrangenum);

clear psalnum tmpnum tempnum presnum i ii ind j k outofrangenum NaNnum ans

badnum = badnum(6:end)

fnumb = []; timeb = [];

badnum(i)
ind = find(fnum == badnum(i));
figure(1); clf; plot(temp2(:,ind),-pres2(:,ind),'.-'); grid on;
figure(2); clf; plot(psal2(:,ind),-pres2(:,ind),'.-'); grid on;
tmpt = temp2(:,ind); tmps = psal2(:,ind); tmpp = pres2(:,ind); tmpj = datestr(time(ind));
if sum(badnum(i) == grayfloat) > 0
  disp('greylisted');
  ii = find(badnum(i) == grayfloat);
  disp(param(ii)); disp(errtype(ii)); disp(start_date(ii)); disp(end_date(ii));
end

i = 1;
% set bad values to Nan
temp2(1,ind(24)) = NaN; 
psal2(isnan(temp2)) = NaN;
fnumb = [fnumb; fnum(ind(24))];
timeb = [timeb; time(ind(24))];

i = 2;
% set bad values to Nan
temp2(1,ind(26)) = NaN; 
temp2(2:3,ind(35)) = NaN;
psal2(isnan(temp2)) = NaN;
fnumb = [fnumb; fnum(ind(26))];
timeb = [timeb; time(ind(26))];

num_lowS = [fnumb timeb];
save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_lowS -APPEND

clear num_lowS fnumb timeb indb ind ans i ii tmpj tmpp tmps tmpt

% for ii = 1:length(ind)
%     ii
%     figure(1); clf; plot(psal2(ind,:)',-pres2(ind,:)','-','color',[.9 .9 .9]);
%     hold on; plot(psal2(ind(ii),:)',-pres2(ind(ii),:)','r.-'); grid on;
%     figure(2); clf;
%     axesm('mercator','MapLatLimit',[min(latitude(ind))-5 max(latitude(ind))+5],'MapLonLimit',[min(longitude(ind))-5 max(longitude(ind))+5])
%     plotm(latitude(ind),longitude(ind),'.-'); hold on;
%     plotm(latitude(ind(ii)),longitude(ind(ii)),'r.');
%     geoshow('landareas.shp','FaceColor',[.9 .9 .9])
%     setm(gca, 'grid','on','meridianlabel','on','parallellabel','on','plinelocation',10,'mlinelocation',10,'plabellocation',10,'mlabellocation',10)
%     axis tight;
%     figure(3); clf; plot(temp2(ind,:)',-pres2(ind,:)','-','color',[.9 .9 .9]);
%     hold on; plot(temp2(ind(ii),:)',-pres2(ind(ii),:)','r.-'); grid on;
%     pause
% end

% delete profiles that are all NaNs
clear ind indb; k = 1; l = 1;
for i = 1:length(fnum)
    t = temp2(:,i); s = psal2(:,i); p = pres2(:,i);
    ss = sum(~isnan(s)); tt = sum(~isnan(t)); pp = sum(~isnan(p));
    
    if ~(tt == 0 || ss == 0 || pp == 0)
        ind(k) = i; k = k+1; %#ok<SAGROW>
    else
        indb(l) = i; l = l+1; %#ok<SAGROW>
    end
end
clear t s p ss tt pp k i l
ind = ind'; indb = indb';

temp2 = temp2(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres2 = pres2(:,ind);
psal2 = psal2(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
pqc = pqc(:,ind); pres = pres(:,ind); psal = psal(:,ind); sqc = sqc(:,ind); temp = temp(:,ind); tqc = tqc(:,ind);
pres_adj = pres_adj(:,ind); pres_err = pres_err(:,ind); pres_qc = pres_qc(:,ind);
psal_adj = psal_adj(:,ind); psal_err = psal_err(:,ind); psal_qc = psal_qc(:,ind);
temp_adj = temp_adj(:,ind); temp_err = temp_err(:,ind); temp_qc = temp_qc(:,ind);
datamode = datamode(ind); time_adj = time_adj(ind);
clear ind indb 

%% find profiles with repeat pressures
clear ind1; k = 1;
for i = 1:length(fnum)
    p = pres2(:,i);
    for j = 1:length(p)
        if ~isnan(p(j))
            if sum(p == p(j)) > 1
                ind1(k) = i; %#ok<SAGROW>
                k = k+1;
            end
        end
    end
end
clear p i j k
ind1 = ind1';

% individual float numbers of profiles with repeat pressures
k = 1; ind = 0;
for i = 1:length(ind1)
    if fnum(ind1(i)) ~= ind
        ind(k) = fnum(ind1(i)); %#ok<SAGROW>
        k = k+1;
    end
end
ind = ind'; clear i k
length(ind)

% for i = 1:length(ind)
%     ii = find(fnum == ind(i));
%     plot(temp2(:,ii),-pres2(:,ii),'.-');
%     ind(i)
%     pres2(1:2,ii)
%     pause
% end
% clear ii i

i = 1;

i = i+1;
ind(i)
ii = find(fnum == ind(i));
figure(1); clf; plot(temp2(:,ii),-pres2(:,ii),'.-'); grid on;
figure(2); clf; plot(psal2(:,ii),-pres2(:,ii),'.-'); grid on;
tmpt = temp2(:,ii); tmps = psal2(:,ii); tmpp = pres2(:,ii); tmpj = datestr(time(ii));
if sum(ind(i) == grayfloat) > 0
  disp('greylisted');
  jj = find(ind(i) == grayfloat);
  disp(param(jj)); disp(errtype(jj)); disp(start_date(jj)); disp(end_date(jj));
end

clear ind2 ind3;
k = 1;
for iii = 1:length(ii)
p = pres2(:,ii(iii));
    for j = 1:length(p)
        if ~isnan(p(j))
            if sum(p == p(j)) > 1
                ind2(k) = iii;  %#ok<SAGROW>
                ind3(k) = j; %#ok<SAGROW>
                k = k+1;
            end
        end
    end
end
disp(ind2)

temp2(ind3(1),ii(ind2(1))) = NaN;

ind = find(isnan(pres2));
temp2(ind) = NaN;
psal2(ind) = NaN;

ind = find(isnan(psal2));
temp2(ind) = NaN;
pres2(ind) = NaN;

ind = find(isnan(temp2));
pres2(ind) = NaN;
psal2(ind) = NaN;

% find profiles with increasing pressures 
k = 1; ind1 = [];
for i = 1:length(fnum)
    p = pres2(:,i);
    p = p(~isnan(p));
    pd = p(2:end)-p(1:end-1);
    if sum(pd < 0) > 0
        ind1(k) = i; %#ok<SAGROW>
        k = k+1;
    end
end
clear ind i k p pd s t j
ind1 = ind1';

% individual float numbers of profiles with increasing pressures
k = 1; ind = 0;
for i = 1:length(ind1)
    if fnum(ind1(i)) ~= ind
        ind(k) = fnum(ind1(i)); %#ok<SAGROW>
        k = k+1;
    end
end
ind = ind'; clear i k
length(ind)

for i = 1:length(ind)
    ii = find(fnum == ind(i));
    plot(temp2(:,ii),-pres2(:,ii),'.-');
    ind(i)
    pres2(1:2,ii)
    pause
end
clear ii i

i = 1;
pres2(12,ii(2)) = NaN; 
pres2(39,ii(32)) = NaN; 
pres2(31,ii(38)) = NaN; 

ind(i)
ii = find(fnum == ind(i));
figure(1); clf; plot(temp2(:,ii),-pres2(:,ii),'.-'); grid on;
figure(2); clf; plot(psal2(:,ii),-pres2(:,ii),'.-'); grid on;
tmpt = temp2(:,ii); tmps = psal2(:,ii); tmpp = pres2(:,ii); tmpj = datestr(time(ii));
if sum(ind(i) == grayfloat) > 0
  disp('greylisted');
  jj = find(ind(i) == grayfloat);
  disp(param(jj)); disp(errtype(jj)); disp(start_date(jj)); disp(end_date(jj));
end

k = 1; clear ind2
for iii = 1:length(ii)
    p = pres2(:,ii(iii));
    p = p(~isnan(p));
    pd = p(2:end)-p(1:end-1);
    if sum(pd < 0) > 0
        ind2(k) = iii; %#ok<SAGROW>
        k = k+1;
    end
end
disp(ind2)

ind = find(isnan(pres2));
temp2(ind) = NaN;
psal2(ind) = NaN;

ind = find(isnan(psal2));
temp2(ind) = NaN;
pres2(ind) = NaN;

ind = find(isnan(temp2));
pres2(ind) = NaN;
psal2(ind) = NaN;

clear ans i ii ind ind1 jj badnum tmpp tmpt tmpj tmps
clear iii ind2 ind3

%% plot all profiles and look for errant ones
figure(1); clf
plot(temp2,-pres2,'.-'); grid on; axis tight;

figure(2); clf 
plot(psal2,-pres2,'.-'); grid on; axis tight;

[i,j] = find(pres2 > 950 & psal2 < 34.05 & pres2 < 1250);
for k = 1:length(i)
    pres2(i(k),j(k)) = NaN;
end

[i,j] = find(pres2 > 1800 & psal2 < 34.35);
for k = 1:length(i)
    pres2(i(k),j(k)) = NaN;
end


% i = 6;
% temp2(:,ind(12)) = NaN;
% psal2(isnan(temp2)) = NaN;
% fnumb = [fnumb; fnum(ind(12))]; timeb = [timeb; time(ind(12))];
% 
% badnum(i)
% ind = find(fnum == badnum(i));
% figure(1); clf; plot(temp2(:,ind),-pres2(:,ind),'.-'); grid on;
% figure(2); clf; plot(psal2(:,ind),-pres2(:,ind),'.-'); grid on;
% tmpt = temp2(:,ind); tmps = psal2(:,ind); tmpp = pres2(:,ind); tmpj = datestr(time(ind));
% if sum(badnum(i) == grayfloat) > 0
%   disp('greylisted');
%   ii = find(badnum(i) == grayfloat);
%   disp(param(ii)); disp(errtype(ii)); disp(start_date(ii)); disp(end_date(ii));
% end
% 
% num_errant = [fnumb timeb];
% save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_errant -APPEND
% clear num_errant fnumb timeb i ii ind indnum j jj badnum ans tmpj tmpp tmps tmpt


% need t, s, and p at every point - make NaN's consistent across all data
ind = find(isnan(pres2));
temp2(ind) = NaN;
psal2(ind) = NaN;

ind = find(isnan(psal2));
temp2(ind) = NaN;
pres2(ind) = NaN;

ind = find(isnan(temp2));
pres2(ind) = NaN;
psal2(ind) = NaN;

clear ind

% delete profiles that are all NaNs
clear ind indb; k = 1; l = 1;
for i = 1:length(fnum)
    t = temp2(:,i); s = psal2(:,i); p = pres2(:,i);
    ss = sum(~isnan(s)); tt = sum(~isnan(t)); pp = sum(~isnan(p));
    
    if ~(tt == 0 || ss == 0 || pp == 0)
        ind(k) = i; k = k+1; %#ok<SAGROW>
    else
        indb(l) = i; l = l+1; %#ok<SAGROW>
    end
end
clear t s p ss tt pp k i l
ind = ind'; indb = indb';

temp2 = temp2(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres2 = pres2(:,ind);
psal2 = psal2(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
pqc = pqc(:,ind); pres = pres(:,ind); psal = psal(:,ind); sqc = sqc(:,ind); temp = temp(:,ind); tqc = tqc(:,ind);
pres_adj = pres_adj(:,ind); pres_err = pres_err(:,ind); pres_qc = pres_qc(:,ind);
psal_adj = psal_adj(:,ind); psal_err = psal_err(:,ind); psal_qc = psal_qc(:,ind);
temp_adj = temp_adj(:,ind); temp_err = temp_err(:,ind); temp_qc = temp_qc(:,ind);
datamode = datamode(ind); time_adj = time_adj(ind);
clear ind indb

clear ans badnum end_date errtype grayfloat ii ind ind1 j jj kk num_* param
clear tmpj tmpp tmps tmpt start_date
% 
%% maximum pressure
maxp = (max(abs(pres2),[],1))';

% find & remove profiles with maxp < 800 db
indp = find(maxp < 800);
per800 = 100*length(indp)/length(maxp)

num_800_3 = [fnum(indp) time(indp)];

ind = find(maxp >= 800);
temp2 = temp2(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres2 = pres2(:,ind);
psal2 = psal2(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
pqc = pqc(:,ind); pres = pres(:,ind); psal = psal(:,ind); sqc = sqc(:,ind); temp = temp(:,ind); tqc = tqc(:,ind);
pres_adj = pres_adj(:,ind); pres_err = pres_err(:,ind); pres_qc = pres_qc(:,ind);
psal_adj = psal_adj(:,ind); psal_err = psal_err(:,ind); psal_qc = psal_qc(:,ind);
temp_adj = temp_adj(:,ind); temp_err = temp_err(:,ind); temp_qc = temp_qc(:,ind);
datamode = datamode(ind); time_adj = time_adj(ind);

clear ans ind indp per800 t2 

% % minimum pressure
% minp = (min(abs(pres2),[],1))';
% 
% % find & remove profiles with maxp < 800 db
% indp = find(minp > 400);
% per400 = 100*length(indp)/length(minp)
% 
% num_400 = [fnum(indp) time(indp)];
% 
% ind = find(minp <= 400);
% temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
% psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
% temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
% position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
% prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
% temp2 = temp2(:,ind); psal2 = psal2(:,ind); pres2 = pres2(:,ind);
% tqc = tqc(:,ind); pqc = pqc(:,ind); sqc = sqc(:,ind);
% 
% clear ans ind indp per800 t2 minp maxp per400

save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_800_3 -APPEND
clear num_800_3 num_400 maxp

clear ind3 jjj outofrangenum psalnum presnum tempnum tmpnum

%% remove fsi sensors
ind = find((type == 842 | (type == 847 | (type == 852 | type == 857))));

k = 1; indnum = 0;
for i = 1:length(fnum(ind))
    if fnum(ind(i)) ~= indnum
        indnum(k) = fnum(ind(i)); %#ok<SAGROW>
        k = k+1;
    end
end
indnum = indnum';
clear i k

cd /home/alison/Matlab/M' files'/Global' Mapping 2'/

for i = 1:length(indnum)
    fn = num2str(indnum(i));
    
    filename = strcat('/home/alison/Data/Argo Profiles/Float/North Pacific/','Profiles/',fn,'_prof.nc');
    nc = netcdf.open(char(filename),'NC_NOWRITE');
    pi = ncchar(nc,'PI_NAME');
    PI(i,:) = (pi(:,1))';
    netcdf.close(nc);
end

indnum2 = NaN(length(indnum),1);
for i = 1:length(indnum)
   
    if strcmp(PI(i,1:5),'BRECK')
        indnum2(i) = 1;
    else indnum2(i) = 0;
    end
end

indnum = indnum(indnum2==0);

for i = 1:length(indnum)

    ind = find(fnum == indnum(i));
    time(ind) = NaN;
end

ind = find(~isnan(time));

temp2 = temp2(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres2 = pres2(:,ind);
psal2 = psal2(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
pqc = pqc(:,ind); pres = pres(:,ind); psal = psal(:,ind); sqc = sqc(:,ind); temp = temp(:,ind); tqc = tqc(:,ind);
pres_adj = pres_adj(:,ind); pres_err = pres_err(:,ind); pres_qc = pres_qc(:,ind);
psal_adj = psal_adj(:,ind); psal_err = psal_err(:,ind); psal_qc = psal_qc(:,ind);
temp_adj = temp_adj(:,ind); temp_err = temp_err(:,ind); temp_qc = temp_qc(:,ind);
datamode = datamode(ind); time_adj = time_adj(ind);

clear PI ans filename fn i ind indnum indnum2 nc pi j


%% %%%%%%%%%%%%%%%
save('/home/alison/Matlab/Mat files/Global Mapping 2/North Pacific/prof_2008_qc_all.mat');   

clear juld_qc percentNaN position_qc pqc pres pres_qc prof_* psal psal_qc sqc dir_name indnum temp temp_qc tqc
pres = pres2; psal = psal2; temp = temp2; clear pres2 psal2 temp2 pres_adj pres_qc_adj psal_adj psal_qc_adj temp_adj temp_qc_adj time_adj

save('/home/alison/Matlab/Mat files/Global Mapping 2/North Pacific/prof_2008_qc.mat');

clear; load /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat
save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats_2008.mat

%% quality control

% 	
clear; pack
load /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/prof_2009.mat;
load /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/prof_adj_2009.mat;
load /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/prof_err_2009.mat;

pres_qc = pres_qc';
psal_qc = psal_qc';
temp_qc = temp_qc';
pres_qc_adj = pres_qc_adj';
psal_qc_adj = psal_qc_adj';
temp_qc_adj = temp_qc_adj';

longitude(longitude < 0) = longitude(longitude < 0) + 360;

% plot locations of data
figure(1); clf
axesm('mercator','MapLatLimit',[-5 80],'MapLonLimit',[100 280])
plotm(latitude,longitude,'.');
geoshow('landareas.shp','FaceColor',[.9 .9 .9])
setm(gca, 'grid','on','meridianlabel','on','parallellabel','on','plinelocation',5,'mlinelocation',5,'plabellocation',5,'mlabellocation',5)
axis tight

figure(2); clf
plot(longitude,latitude,'.'); grid on;

% lat/lon range
ind = find(longitude <= 283 & longitude >= 120.5 & latitude <= 63 & latitude >= 0);
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

% remove sea of japan profiles
ind = find(~((longitude < 138 & latitude < 44.1 & latitude > 34.5) | (longitude < 141.1 & latitude < 45.1 & latitude > 37) | (longitude < 141.7 & latitude < 45.3 & latitude > 44))); 
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

ind = find(~(longitude > 250 & latitude > 33.5 & latitude < 57));
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

ind = find(~(longitude < 121 & latitude > 18 & latitude < 18.5));
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

ind = find(~(longitude < 125 & latitude > 11 & latitude < 11.3));
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

ind = find(~(longitude < 125.8 & latitude > 2 & latitude < 6.4));
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

ind = find(~(longitude > 243.64 & latitude > 31.21 & latitude < 31.235));
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

% plot locations of data
figure(1); clf
axesm('mercator','MapLatLimit',[-5 80],'MapLonLimit',[100 280])
plotm(latitude,longitude,'.');
geoshow('landareas.shp','FaceColor',[.9 .9 .9])
setm(gca, 'grid','on','meridianlabel','on','parallellabel','on','plinelocation',10,'mlinelocation',10,'plabellocation',10,'mlabellocation',10)
axis tight

% look at maximum pressures
maxp = (max(abs(pres),[],1))';

% find & remove profiles with maxp < 800 db
indp = find(maxp < 800);
per800 = 100*length(indp)/length(maxp) %#ok<*NOPTS,*NASGU>

num_800 = [fnum(indp) time(indp)];

ind = find(maxp >= 800);
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

clear ans ind indp per800 t2 

save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_800
clear num_800 maxp

% delete all Nan profiles
clear ind; k = 1; 
for i = 1:length(fnum)
    t = temp(:,i); s = psal(:,i); p = pres(:,i);
    ss = sum(~isnan(s)); tt = sum(~isnan(t)); pp = sum(~isnan(p));
    
    if ~(tt == 0 || ss == 0 || pp == 0)
        ind(k) = i; k = k+1; %#ok<NAGROW>
    end
end
clear t s p ss tt pp k i
ind = ind';
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

indt = find(sum(~isnan(temp),2) == 0); %#ok<MXFND>
indp = find(sum(~isnan(pres),2) == 0); %#ok<MXFND>
inds = find(sum(~isnan(psal),2) == 0); %#ok<MXFND>
ind = max([min(indt) min(indp) min(inds)]);
temp = temp(1:ind,:); psal = psal(1:ind,:); pres = pres(1:ind,:);
temp_qc = temp_qc(1:ind,:); pres_qc = pres_qc(1:ind,:); psal_qc = psal_qc(1:ind,:);
temp_adj = temp_adj(1:ind,:); psal_adj = psal_adj(1:ind,:); pres_adj = pres_adj(1:ind,:);
temp_qc_adj = temp_qc_adj(1:ind,:); pres_qc_adj = pres_qc_adj(1:ind,:); psal_qc_adj = psal_qc_adj(1:ind,:);
temp_err = temp_err(1:ind,:); psal_err = psal_err(1:ind,:); pres_err = pres_err(1:ind,:);
clear indt indp inds ind

save /home/alison/Matlab/M'at files'/Global' Mapping 2'/North' Pacific'/tmp.mat

% graylist floats
cd('/home/alison/Matlab/M files/Global Mapping 2/');
[grayfloat,~,start_date,end_date,errtype] = graylist('/home/alison/Matlab/Mat files/Global Mapping 2/North Pacific/');


%% %%%% quality control flags %%%%

% examine flags on temperature data
[mm,nn] = size(temp);
tqc = NaN(size(temp));
for j = 1:nn
    for i = 1:mm
        if ~isnan(temp(i,j))
            if (strcmp(temp_qc_adj(i,j),'1') || strcmp(temp_qc_adj(i,j),'0') || strcmp(temp_qc_adj(i,j),'2') || strcmp(temp_qc_adj(i,j),' ') || double(temp_qc_adj(i,j))==0)
                tqc(i,j) = 1;
            else
                tqc(i,j) = 0;
            end
        end
    end
end

% examine flags on salinity data
[mm,nn] = size(psal);
sqc = NaN(size(psal));
for j = 1:nn
    for i = 1:mm
        if ~isnan(psal(i,j))
            if (strcmp(psal_qc_adj(i,j),'1') || strcmp(psal_qc_adj(i,j),'0') || strcmp(psal_qc_adj(i,j),'2') || strcmp(psal_qc_adj(i,j),' ') || double(psal_qc_adj(i,j))==0)
                sqc(i,j) = 1;
            else
                sqc(i,j) = 0;
            end
        end
    end
end

% examine flags on pressure data
[mm,nn] = size(pres);
pqc = NaN(size(pres));
for j = 1:nn
    for i = 1:mm
        if ~isnan(pres(i,j))
            if (strcmp(pres_qc_adj(i,j),'1') || strcmp(pres_qc_adj(i,j),'0') || strcmp(pres_qc_adj(i,j),'2') || strcmp(pres_qc_adj(i,j),' ') || double(pres_qc_adj(i,j))==0)
                pqc(i,j) = 1;
            else
                pqc(i,j) = 0;
            end
        end
    end
end

clear mm nn i j


% display number of profiles with one or more bad qc flag
num_bad_temp = length(find(sum(tqc == 0,1) > 0))
num_bad_pres = length(find(sum(pqc == 0,1) > 0))
num_bad_psal = length(find(sum(sqc == 0,1) > 0))
clear num_bad_*

% make 2nd variable set with data with qc ~= 0,1,2,' ' set to NaN
temp2 = temp_adj; temp2(tqc == 0) = NaN;
psal2 = psal_adj; psal2(sqc == 0) = NaN;
pres2 = pres_adj; pres2(pqc == 0) = NaN;

% need t, s, and p at every point - make NaN's consistent across all data
ind = find(isnan(pres2));
temp2(ind) = NaN;
psal2(ind) = NaN;

ind = find(isnan(psal2));
temp2(ind) = NaN;
pres2(ind) = NaN;

ind = find(isnan(temp2));
pres2(ind) = NaN;
psal2(ind) = NaN;

% delete profiles that are all NaNs
clear ind indb; k = 1; l = 1;
for i = 1:length(fnum)
    t = temp2(:,i); s = psal2(:,i); p = pres2(:,i);
    ss = sum(~isnan(s)); tt = sum(~isnan(t)); pp = sum(~isnan(p));
    
    if ~(tt == 0 || ss == 0 || pp == 0)
        ind(k) = i; k = k+1; %#ok<SAGROW>
    else
        indb(l) = i; l = l+1; %#ok<SAGROW>
    end
end
clear t s p ss tt pp k i l
ind = ind'; indb = indb';

num_badQC = [fnum(indb) time(indb)];
save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_badQC -APPEND
clear indb num_badQC

temp2 = temp2(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres2 = pres2(:,ind);
psal2 = psal2(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
pqc = pqc(:,ind); pres = pres(:,ind); psal = psal(:,ind); sqc = sqc(:,ind); temp = temp(:,ind); tqc = tqc(:,ind);
pres_adj = pres_adj(:,ind); pres_err = pres_err(:,ind); pres_qc = pres_qc(:,ind);
psal_adj = psal_adj(:,ind); psal_err = psal_err(:,ind); psal_qc = psal_qc(:,ind);
temp_adj = temp_adj(:,ind); temp_err = temp_err(:,ind); temp_qc = temp_qc(:,ind);
datamode = datamode(ind); time_adj = time_adj(ind);

clear ind

% find & remove profiles with maxp < 800 db
maxp = (max(abs(pres2),[],1))';
indp = find(maxp < 800);
per800 = 100*length(indp)/length(maxp) %#ok<*NOPTS,*NASGU>

num_800_2 = [fnum(indp) time(indp)];

ind = find(maxp >= 800);
temp2 = temp2(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres2 = pres2(:,ind);
psal2 = psal2(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
pqc = pqc(:,ind); pres = pres(:,ind); psal = psal(:,ind); sqc = sqc(:,ind); temp = temp(:,ind); tqc = tqc(:,ind);
pres_adj = pres_adj(:,ind); pres_err = pres_err(:,ind); pres_qc = pres_qc(:,ind);
psal_adj = psal_adj(:,ind); psal_err = psal_err(:,ind); psal_qc = psal_qc(:,ind);
temp_adj = temp_adj(:,ind); temp_err = temp_err(:,ind); temp_qc = temp_qc(:,ind);
datamode = datamode(ind); time_adj = time_adj(ind);
clear ind

clear ans ind indp per800 t2 

save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_800_2
clear num_800_2 maxp

% % compute percentage of good values in each profile
% [~,nn] = size(pres2);
% percentNaN = NaN(nn,3);
% for i = 1:nn
%     p = pqc(:,i); t = tqc(:,i); s = sqc(:,i);
%     np = sum(p == 0); nt = sum(t == 0); ns = sum(s == 0); 
%     npp = sum(~isnan(p)); ntt = sum(~isnan(t)); nss = sum(~isnan(s));
%     percentNaN(i,:) = 100*[np./npp nt/ntt ns/nss];
% end
% clear i p np lastp mm nn npp t nt ntt s ns nss
% 
% ind = find(percentNaN(:,1) > 25 | percentNaN(:,2) > 25 | percentNaN(:,3) > 25);
% for ii = 1:length(ind)
%     fnum(ind(ii))
%     figure(1); clf
%     plot(temp(:,ind(ii)),-pres(:,ind(ii)),'m.-'); hold on;
%     plot(temp2(:,ind(ii)),-pres2(:,ind(ii)),'r.-');
%     disp(temp_qc(:,ind(ii))');
%     disp(pres_qc(:,ind(ii))');
%     disp(psal_qc(:,ind(ii))');
%     figure(2); clf
%     plot(psal(:,ind(ii)),-pres(:,ind(ii)),'c.-'); hold on;
%     plot(psal2(:,ind(ii)),-pres2(:,ind(ii)),'b.-');
%     pause
% end
% clear ii ans
% 
% num_percentNaN = [fnum(ind) time(ind)];
% save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_percentNaN -APPEND
% clear num_percentNaN
% 
% temp2(:,ind) = NaN; pres2(:,ind) = NaN; psal2(:,ind) = NaN;

% find & delete short profiles
n = sum(~isnan(temp2),1)';
% for i = 1:10
%     disp(i); ind = find(n == i); length(ind)
%     figure(1); clf
%     plot(temp2(:,ind),-pres2(:,ind),'.-');
%     pause
% end

ind = find(n <= 10 & n > 0); length(ind)
num_shortprof = [fnum(ind) time(ind)];
save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_shortprof -APPEND
clear num_shortprof

temp2(:,ind) = NaN; pres2(:,ind) = NaN; psal2(:,ind) = NaN;
clear n i ind ans

% delete profiles that are all NaNs
clear ind indb; k = 1; l = 1;
for i = 1:length(fnum)
    t = temp2(:,i); s = psal2(:,i); p = pres2(:,i);
    ss = sum(~isnan(s)); tt = sum(~isnan(t)); pp = sum(~isnan(p));
    
    if ~(tt == 0 || ss == 0 || pp == 0)
        ind(k) = i; k = k+1; %#ok<SAGROW>
    else
        indb(l) = i; l = l+1; %#ok<SAGROW>
    end
end
clear t s p ss tt pp k i l
ind = ind'; indb = indb';

temp2 = temp2(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres2 = pres2(:,ind);
psal2 = psal2(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
pqc = pqc(:,ind); pres = pres(:,ind); psal = psal(:,ind); sqc = sqc(:,ind); temp = temp(:,ind); tqc = tqc(:,ind);
pres_adj = pres_adj(:,ind); pres_err = pres_err(:,ind); pres_qc = pres_qc(:,ind);
psal_adj = psal_adj(:,ind); psal_err = psal_err(:,ind); psal_qc = psal_qc(:,ind);
temp_adj = temp_adj(:,ind); temp_err = temp_err(:,ind); temp_qc = temp_qc(:,ind);
datamode = datamode(ind); time_adj = time_adj(ind);
clear ind indb 

save /home/alison/Matlab/M'at files'/Global' Mapping 2'/North' Pacific'/tmp.mat

%% remove greylist floats
k = 1;
badfloats = 0; type_err = cell(1,1); bad_start = 0; bad_end = 0;
for i = 1:length(grayfloat)
    if sum(grayfloat(i) == fnum) > 0
    if grayfloat(i) ~= badfloats
        badfloats(k,1) = grayfloat(i); %#ok<SAGROW>
        type_err(k,1) = errtype(i);
        bad_start(k,1) = start_date(i); %#ok<SAGROW>
        bad_end(k,1) = end_date(i); %#ok<SAGROW>
        k = k+1;
    end
    end
end
clear k i

% ind = find(~strcmp(type_err,'potential problem with bin assignment'));
% badfloats = badfloats(ind);
% type_err = type_err(ind);
% bad_start = bad_start(ind);
% bad_end = bad_end(ind);

badstart = cellstr(num2str(bad_start));
badend = cellstr(num2str(bad_end));

for i = 1:length(badstart)
    bs = char(badstart(i));
    bs = datenum([sscanf(bs(1:4),'%f'),sscanf(bs(5:6),'%f'),sscanf(bs(7:8),'%f')]);
    bad_st(i,1) = bs; %#ok<SAGROW>
    be = char(badend(i));
    if ~strcmp(be,'     NaN') && ~strcmp(be,'NaN')
        be = datenum([sscanf(be(1:4),'%f'),sscanf(be(5:6),'%f'),sscanf(be(7:8),'%f')]);
        bad_end(i) = be;
    else
        bad_end(i) = NaN;
    end
end

fnumb = []; timeb = [];
for j = 1:length(badfloats)
    if isnan(bad_end(j))
        ind = find(~((fnum == badfloats(j))&(time > bad_st(j))));
        indb = find((fnum == badfloats(j)) & (time > bad_st(j)));
    else
        ind = find(~(fnum == badfloats(j) & time > bad_st(j) & time < bad_end(j)));
        indb = find(fnum == badfloats(j) & time > bad_st(j) & time < bad_end(j));
    end
fnumb = [fnumb; fnum(indb)]; %#ok<AGROW>
timeb = [timeb; time(indb)]; %#ok<AGROW>

temp2 = temp2(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres2 = pres2(:,ind);
psal2 = psal2(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
pqc = pqc(:,ind); pres = pres(:,ind); psal = psal(:,ind); sqc = sqc(:,ind); temp = temp(:,ind); tqc = tqc(:,ind);
pres_adj = pres_adj(:,ind); pres_err = pres_err(:,ind); pres_qc = pres_qc(:,ind);
psal_adj = psal_adj(:,ind); psal_err = psal_err(:,ind); psal_qc = psal_qc(:,ind);
temp_adj = temp_adj(:,ind); temp_err = temp_err(:,ind); temp_qc = temp_qc(:,ind);
datamode = datamode(ind); time_adj = time_adj(ind);

end

num_gray = [fnumb timeb];
save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_gray -APPEND

clear num_gray fnumb timeb indb ind
clear bad* end_date j ind param start_date be bs i 
clear ans errtype grayfloat ii percentNaN type_err

% graylist floats
cd('/home/alison/Matlab/M files/Global Mapping 2/');
[grayfloat,param,start_date,end_date,errtype] = graylist('/home/alison/Matlab/Mat files/Global Mapping 2/North Pacific/');

save /home/alison/Matlab/M'at files'/Global' Mapping 2'/North' Pacific'/tmp.mat

%% look for out-of-range pressures
figure(1); clf; imagesc(pres); colorbar;
figure(2); clf; imagesc(pres2); colorbar

outofrangenum = [];

ind = find(pres > 2200); length(ind)
[i,j] = find(pres2 > 2200); length(i)
tempnum = fnum(j);
presnum = 0; k = 1;
for ii = 1:length(tempnum)
if tempnum(ii) ~= presnum
    presnum(k,1) = tempnum(ii); %#ok<SAGROW>
    k = k+1;
end
end
outofrangenum = [outofrangenum; presnum];

% for ii = 1:length(i)
%     figure(1); clf
%     plot(temp2(i(ii),:),-pres2(i(ii),:),'r.-'); hold on
%     plot(psal2(i(ii),:),-pres2(i(ii),:),'b.-');
%     pres2(i(ii),1:30)
%     pause
% end

ind = find(pres < 0); length(ind)
[~,j] = find(pres2 < 0); length(j)
tempnum = fnum(j);
presnum = 0; k = 1;
for ii = 1:length(tempnum)
if tempnum(ii) ~= presnum
    presnum(k,1) = tempnum(ii); %#ok<SAGROW>
    k = k+1;
end
end
outofrangenum = [outofrangenum; presnum];

% for ii = 1:length(i)
%     figure(1); clf
%     plot(temp(i(ii),:),-pres(i(ii),:),'r.-'); hold on
%     plot(psal(i(ii),:),-pres(i(ii),:),'b.-');
%     pres2(i(ii),1:30)
%     pause
% end

% look for out-of-range temperatures
figure(1); clf; imagesc(temp); colorbar;
figure(2); clf; imagesc(temp2); colorbar

ind = find(temp > 40); length(ind)
[~,j] = find(temp2 > 40); length(j)
tempnum = fnum(j);
tmpnum = 0; k = 1;
for ii = 1:length(tempnum)
if tempnum(ii) ~= tmpnum
    tmpnum(k,1) = tempnum(ii); %#ok<SAGROW>
    k = k+1;
end
end
outofrangenum = [outofrangenum; tmpnum];

% for ii = 1:length(i)
%     figure(1); clf
%     plot(temp2(i(ii),:),-pres2(i(ii),:),'r.-'); hold on
%     plot(psal2(i(ii),:),-pres2(i(ii),:),'b.-');
%     temp2(i(ii),1:30)
%     pause
% end

ind = find(temp <= -3); length(ind)
[~,j] = find(temp2 <= -3); length(j)
tempnum = fnum(j);
tmpnum = 0; k = 1;
for ii = 1:length(tempnum)
if tempnum(ii) ~= tmpnum
    tmpnum(k,1) = tempnum(ii); %#ok<SAGROW>
    k = k+1;
end
end
outofrangenum = [outofrangenum; tmpnum];

% for ii = 1:length(i)
%     figure(1); clf
%     plot(temp2(i(ii),:),-pres2(i(ii),:),'r.-'); hold on
%     plot(psal2(i(ii),:),-pres2(i(ii),:),'b.-');
%     temp2(i(ii),1:30)
%     pause
% end

% look for out-of-range salinities
figure(1); clf; imagesc(psal); colorbar;
figure(2); clf; imagesc(psal2); colorbar

ind = find(psal > 39); length(ind)
[~,j] = find(psal2 > 39); length(j)
tempnum = fnum(j);
psalnum = 0; k = 1;
for ii = 1:length(tempnum)
if tempnum(ii) ~= psalnum
    psalnum(k,1) = tempnum(ii); %#ok<SAGROW>
    k = k+1;
end
end
outofrangenum = [outofrangenum; psalnum];

% for ii = 1:length(i)
%     figure(1); clf
% %     plot(temp2(i(ii),:),-pres2(i(ii),:),'r.-'); hold on
%     plot(psal2(:,j(ii)),-pres2(:,j(ii)),'b.-'); grid on
% %     psal2(i(ii),1:30)
%     pause
% end

ind = find(psal <= 20); length(ind)
[~,j] = find(psal2 <= 20); length(j)
tempnum = fnum(j);
psalnum = 0; k = 1;
for ii = 1:length(tempnum)
if tempnum(ii) ~= psalnum
    psalnum(k,1) = tempnum(ii); %#ok<SAGROW>
    k = k+1;
end
end
outofrangenum = [outofrangenum; psalnum];

% for ii = 1:length(i)
%     figure(1); clf
% %     plot(temp2(i(ii),:),-pres2(i(ii),:),'r.-'); hold on
%     plot(psal2(:,j(ii)),-pres2(:,j(ii)),'b.-'); grid on
% %     psal2(i(ii),1:30)
%     pause
% end

badnum = sort(outofrangenum);

clear psalnum tmpnum tempnum presnum i ii ind j k outofrangenum NaNnum ans

badnum = badnum(5:end)

fnumb = []; timeb = [];

badnum(i)
ind = find(fnum == badnum(i));
figure(1); clf; plot(temp2(:,ind),-pres2(:,ind),'.-'); grid on;
figure(2); clf; plot(psal2(:,ind),-pres2(:,ind),'.-'); grid on;
tmpt = temp2(:,ind); tmps = psal2(:,ind); tmpp = pres2(:,ind); tmpj = datestr(time(ind));
if sum(badnum(i) == grayfloat) > 0
  disp('greylisted');
  ii = find(badnum(i) == grayfloat);
  disp(param(ii)); disp(errtype(ii)); disp(start_date(ii)); disp(end_date(ii));
end

i = 1;
% set bad values to Nan
[ii,j] = find(psal2(:,ind) < 34 & pres2(:,ind) > 600);
for k = 1:length(ii)
    pres2(ii(k),ind(j(k))) = NaN;
end
[ii,j] = find(psal2(:,ind) < 34.5 & pres2(:,ind) > 1850);
for k = 1:length(ii)
    pres2(ii(k),ind(j(k))) = NaN;
end
[ii,j] = find(psal2(:,ind) > 34.8 & pres2(:,ind) > 1200);
for k = 1:length(ii)
    pres2(ii(k),ind(j(k))) = NaN;
end
[ii,j] = find(psal2(:,ind) > 34.3 & pres2(:,ind) > 700 & pres2(:,ind) < 780);
for k = 1:length(ii)
    pres2(ii(k),ind(j(k))) = NaN;
end
[ii,j] = find(psal2(:,ind) < 34.2 & pres2(:,ind) > 920 & pres2(:,ind) < 940);
for k = 1:length(ii)
    pres2(ii(k),ind(j(k))) = NaN;
end
clear ii j k

i = 2;
% set bad values to Nan
temp2(:,ind(29)) = NaN; 
psal2(isnan(temp2)) = NaN;
fnumb = [fnumb; fnum(ind(29))];
timeb = [timeb; time(ind(29))];

i = 3;
% set bad values to Nan
temp2(66:67,ind(19)) = NaN; 
psal2(isnan(temp2)) = NaN;
fnumb = [fnumb; fnum(ind(19))];
timeb = [timeb; time(ind(19))];

num_lowS = [fnumb timeb];
save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_lowS -APPEND

clear num_lowS fnumb timeb indb ind ans i ii tmpj tmpp tmps tmpt

% for ii = 1:length(ind)
%     ii
%     figure(1); clf; plot(psal2(ind,:)',-pres2(ind,:)','-','color',[.9 .9 .9]);
%     hold on; plot(psal2(ind(ii),:)',-pres2(ind(ii),:)','r.-'); grid on;
%     figure(2); clf;
%     axesm('mercator','MapLatLimit',[min(latitude(ind))-5 max(latitude(ind))+5],'MapLonLimit',[min(longitude(ind))-5 max(longitude(ind))+5])
%     plotm(latitude(ind),longitude(ind),'.-'); hold on;
%     plotm(latitude(ind(ii)),longitude(ind(ii)),'r.');
%     geoshow('landareas.shp','FaceColor',[.9 .9 .9])
%     setm(gca, 'grid','on','meridianlabel','on','parallellabel','on','plinelocation',10,'mlinelocation',10,'plabellocation',10,'mlabellocation',10)
%     axis tight;
%     figure(3); clf; plot(temp2(ind,:)',-pres2(ind,:)','-','color',[.9 .9 .9]);
%     hold on; plot(temp2(ind(ii),:)',-pres2(ind(ii),:)','r.-'); grid on;
%     pause
% end

% delete profiles that are all NaNs
clear ind indb; k = 1; l = 1;
for i = 1:length(fnum)
    t = temp2(:,i); s = psal2(:,i); p = pres2(:,i);
    ss = sum(~isnan(s)); tt = sum(~isnan(t)); pp = sum(~isnan(p));
    
    if ~(tt == 0 || ss == 0 || pp == 0)
        ind(k) = i; k = k+1; %#ok<SAGROW>
    else
        indb(l) = i; l = l+1; %#ok<SAGROW>
    end
end
clear t s p ss tt pp k i l
ind = ind'; indb = indb';

temp2 = temp2(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres2 = pres2(:,ind);
psal2 = psal2(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
pqc = pqc(:,ind); pres = pres(:,ind); psal = psal(:,ind); sqc = sqc(:,ind); temp = temp(:,ind); tqc = tqc(:,ind);
pres_adj = pres_adj(:,ind); pres_err = pres_err(:,ind); pres_qc = pres_qc(:,ind);
psal_adj = psal_adj(:,ind); psal_err = psal_err(:,ind); psal_qc = psal_qc(:,ind);
temp_adj = temp_adj(:,ind); temp_err = temp_err(:,ind); temp_qc = temp_qc(:,ind);
datamode = datamode(ind); time_adj = time_adj(ind);
clear ind indb 

%% find profiles with repeat pressures
clear ind1; k = 1;
for i = 1:length(fnum)
    p = pres2(:,i);
    for j = 1:length(p)
        if ~isnan(p(j))
            if sum(p == p(j)) > 1
                ind1(k) = i; %#ok<SAGROW>
                k = k+1;
            end
        end
    end
end
clear p i j k
ind1 = ind1';

% individual float numbers of profiles with repeat pressures
k = 1; ind = 0;
for i = 1:length(ind1)
    if fnum(ind1(i)) ~= ind
        ind(k) = fnum(ind1(i)); %#ok<SAGROW>
        k = k+1;
    end
end
ind = ind'; clear i k
length(ind)

% for i = 1:length(ind)
%     ii = find(fnum == ind(i));
%     plot(temp2(:,ii),-pres2(:,ii),'.-');
%     ind(i)
%     pres2(1:2,ii)
%     pause
% end
% clear ii i

i = 1;

i = i+1;
ind(i)
ii = find(fnum == ind(i));
figure(1); clf; plot(temp2(:,ii),-pres2(:,ii),'.-'); grid on;
figure(2); clf; plot(psal2(:,ii),-pres2(:,ii),'.-'); grid on;
tmpt = temp2(:,ii); tmps = psal2(:,ii); tmpp = pres2(:,ii); tmpj = datestr(time(ii));
if sum(ind(i) == grayfloat) > 0
  disp('greylisted');
  jj = find(ind(i) == grayfloat);
  disp(param(jj)); disp(errtype(jj)); disp(start_date(jj)); disp(end_date(jj));
end

clear ind2 ind3;
k = 1;
for iii = 1:length(ii)
p = pres2(:,ii(iii));
    for j = 1:length(p)
        if ~isnan(p(j))
            if sum(p == p(j)) > 1
                ind2(k) = iii;  %#ok<SAGROW>
                ind3(k) = j; %#ok<SAGROW>
                k = k+1;
            end
        end
    end
end
disp(ind2)

temp2(ind3(1),ii(ind2(1))) = NaN;

ind = find(isnan(pres2));
temp2(ind) = NaN;
psal2(ind) = NaN;

ind = find(isnan(psal2));
temp2(ind) = NaN;
pres2(ind) = NaN;

ind = find(isnan(temp2));
pres2(ind) = NaN;
psal2(ind) = NaN;

% find profiles with increasing pressures 
k = 1; ind1 = [];
for i = 1:length(fnum)
    p = pres2(:,i);
    p = p(~isnan(p));
    pd = p(2:end)-p(1:end-1);
    if sum(pd < 0) > 0
        ind1(k) = i; %#ok<SAGROW>
        k = k+1;
    end
end
clear ind i k p pd s t j
ind1 = ind1';

% individual float numbers of profiles with increasing pressures
k = 1; ind = 0;
for i = 1:length(ind1)
    if fnum(ind1(i)) ~= ind
        ind(k) = fnum(ind1(i)); %#ok<SAGROW>
        k = k+1;
    end
end
ind = ind'; clear i k
length(ind)

for i = 1:length(ind)
    ii = find(fnum == ind(i));
    plot(temp2(:,ii),-pres2(:,ii),'.-');
    ind(i)
    pres2(1:2,ii)
    pause
end
clear ii i

i = 1;
pres2(35,ii(3)) = NaN; 

ind(i)
ii = find(fnum == ind(i));
figure(1); clf; plot(temp2(:,ii),-pres2(:,ii),'.-'); grid on;
figure(2); clf; plot(psal2(:,ii),-pres2(:,ii),'.-'); grid on;
tmpt = temp2(:,ii); tmps = psal2(:,ii); tmpp = pres2(:,ii); tmpj = datestr(time(ii));
if sum(ind(i) == grayfloat) > 0
  disp('greylisted');
  jj = find(ind(i) == grayfloat);
  disp(param(jj)); disp(errtype(jj)); disp(start_date(jj)); disp(end_date(jj));
end

k = 1; clear ind2
for iii = 1:length(ii)
    p = pres2(:,ii(iii));
    p = p(~isnan(p));
    pd = p(2:end)-p(1:end-1);
    if sum(pd < 0) > 0
        ind2(k) = iii; %#ok<SAGROW>
        k = k+1;
    end
end
disp(ind2)

ind = find(isnan(pres2));
temp2(ind) = NaN;
psal2(ind) = NaN;

ind = find(isnan(psal2));
temp2(ind) = NaN;
pres2(ind) = NaN;

ind = find(isnan(temp2));
pres2(ind) = NaN;
psal2(ind) = NaN;

clear ans i ii ind ind1 jj badnum tmpp tmpt tmpj tmps
clear iii ind2 ind3

%% plot all profiles and look for errant ones
figure(1); clf
plot(temp2,-pres2,'.-'); grid on; axis tight;

figure(2); clf 
plot(psal2,-pres2,'.-'); grid on; axis tight;

% [i,j] = find(pres2 > 950 & psal2 < 34.05 & pres2 < 1250);
% for k = 1:length(i)
%     pres2(i(k),j(k)) = NaN;
% end
% 
% [i,j] = find(pres2 > 1800 & psal2 < 34.35);
% for k = 1:length(i)
%     pres2(i(k),j(k)) = NaN;
% end

% i = 6;
% temp2(:,ind(12)) = NaN;
% psal2(isnan(temp2)) = NaN;
% fnumb = [fnumb; fnum(ind(12))]; timeb = [timeb; time(ind(12))];
% 
% badnum(i)
% ind = find(fnum == badnum(i));
% figure(1); clf; plot(temp2(:,ind),-pres2(:,ind),'.-'); grid on;
% figure(2); clf; plot(psal2(:,ind),-pres2(:,ind),'.-'); grid on;
% tmpt = temp2(:,ind); tmps = psal2(:,ind); tmpp = pres2(:,ind); tmpj = datestr(time(ind));
% if sum(badnum(i) == grayfloat) > 0
%   disp('greylisted');
%   ii = find(badnum(i) == grayfloat);
%   disp(param(ii)); disp(errtype(ii)); disp(start_date(ii)); disp(end_date(ii));
% end
% 
% num_errant = [fnumb timeb];
% save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_errant -APPEND
% clear num_errant fnumb timeb i ii ind indnum j jj badnum ans tmpj tmpp tmps tmpt


% need t, s, and p at every point - make NaN's consistent across all data
ind = find(isnan(pres2));
temp2(ind) = NaN;
psal2(ind) = NaN;

ind = find(isnan(psal2));
temp2(ind) = NaN;
pres2(ind) = NaN;

ind = find(isnan(temp2));
pres2(ind) = NaN;
psal2(ind) = NaN;

clear ind

% delete profiles that are all NaNs
clear ind indb; k = 1; l = 1;
for i = 1:length(fnum)
    t = temp2(:,i); s = psal2(:,i); p = pres2(:,i);
    ss = sum(~isnan(s)); tt = sum(~isnan(t)); pp = sum(~isnan(p));
    
    if ~(tt == 0 || ss == 0 || pp == 0)
        ind(k) = i; k = k+1; %#ok<SAGROW>
    else
        indb(l) = i; l = l+1; %#ok<SAGROW>
    end
end
clear t s p ss tt pp k i l
ind = ind'; indb = indb';

temp2 = temp2(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres2 = pres2(:,ind);
psal2 = psal2(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
pqc = pqc(:,ind); pres = pres(:,ind); psal = psal(:,ind); sqc = sqc(:,ind); temp = temp(:,ind); tqc = tqc(:,ind);
pres_adj = pres_adj(:,ind); pres_err = pres_err(:,ind); pres_qc = pres_qc(:,ind);
psal_adj = psal_adj(:,ind); psal_err = psal_err(:,ind); psal_qc = psal_qc(:,ind);
temp_adj = temp_adj(:,ind); temp_err = temp_err(:,ind); temp_qc = temp_qc(:,ind);
datamode = datamode(ind); time_adj = time_adj(ind);
clear ind indb

clear ans badnum end_date errtype grayfloat ii ind ind1 j jj kk num_* param
clear tmpj tmpp tmps tmpt start_date
% 
%% maximum pressure
maxp = (max(abs(pres2),[],1))';

% find & remove profiles with maxp < 800 db
indp = find(maxp < 800);
per800 = 100*length(indp)/length(maxp)

num_800_3 = [fnum(indp) time(indp)];

ind = find(maxp >= 800);
temp2 = temp2(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres2 = pres2(:,ind);
psal2 = psal2(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
pqc = pqc(:,ind); pres = pres(:,ind); psal = psal(:,ind); sqc = sqc(:,ind); temp = temp(:,ind); tqc = tqc(:,ind);
pres_adj = pres_adj(:,ind); pres_err = pres_err(:,ind); pres_qc = pres_qc(:,ind);
psal_adj = psal_adj(:,ind); psal_err = psal_err(:,ind); psal_qc = psal_qc(:,ind);
temp_adj = temp_adj(:,ind); temp_err = temp_err(:,ind); temp_qc = temp_qc(:,ind);
datamode = datamode(ind); time_adj = time_adj(ind);

clear ans ind indp per800 t2 

% % minimum pressure
% minp = (min(abs(pres2),[],1))';
% 
% % find & remove profiles with maxp < 800 db
% indp = find(minp > 400);
% per400 = 100*length(indp)/length(minp)
% 
% num_400 = [fnum(indp) time(indp)];
% 
% ind = find(minp <= 400);
% temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
% psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
% temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
% position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
% prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
% temp2 = temp2(:,ind); psal2 = psal2(:,ind); pres2 = pres2(:,ind);
% tqc = tqc(:,ind); pqc = pqc(:,ind); sqc = sqc(:,ind);
% 
% clear ans ind indp per800 t2 minp maxp per400

save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_800_3 -APPEND
clear num_800_3 num_400 maxp

clear ind3 jjj outofrangenum psalnum presnum tempnum tmpnum

%% remove fsi sensors
ind = find((type == 842 | (type == 847 | (type == 852 | type == 857))));

k = 1; indnum = 0;
for i = 1:length(fnum(ind))
    if fnum(ind(i)) ~= indnum
        indnum(k) = fnum(ind(i)); %#ok<SAGROW>
        k = k+1;
    end
end
indnum = indnum';
clear i k

cd /home/alison/Matlab/M' files'/Global' Mapping 2'/

for i = 1:length(indnum)
    fn = num2str(indnum(i));
    
    filename = strcat('/home/alison/Data/Argo Profiles/Float/North Pacific/','Profiles/',fn,'_prof.nc');
    nc = netcdf.open(char(filename),'NC_NOWRITE');
    pi = ncchar(nc,'PI_NAME');
    PI(i,:) = (pi(:,1))';
    netcdf.close(nc);
end

indnum2 = NaN(length(indnum),1);
for i = 1:length(indnum)
   
    if strcmp(PI(i,1:5),'BRECK')
        indnum2(i) = 1;
    else indnum2(i) = 0;
    end
end

indnum = indnum(indnum2==0);

for i = 1:length(indnum)

    ind = find(fnum == indnum(i));
    time(ind) = NaN;
end

ind = find(~isnan(time));

temp2 = temp2(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres2 = pres2(:,ind);
psal2 = psal2(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
pqc = pqc(:,ind); pres = pres(:,ind); psal = psal(:,ind); sqc = sqc(:,ind); temp = temp(:,ind); tqc = tqc(:,ind);
pres_adj = pres_adj(:,ind); pres_err = pres_err(:,ind); pres_qc = pres_qc(:,ind);
psal_adj = psal_adj(:,ind); psal_err = psal_err(:,ind); psal_qc = psal_qc(:,ind);
temp_adj = temp_adj(:,ind); temp_err = temp_err(:,ind); temp_qc = temp_qc(:,ind);
datamode = datamode(ind); time_adj = time_adj(ind);

clear PI ans filename fn i ind indnum indnum2 nc pi j


%% %%%%%%%%%%%%%%%
save('/home/alison/Matlab/Mat files/Global Mapping 2/North Pacific/prof_2009_qc_all.mat');   

clear juld_qc percentNaN position_qc pqc pres pres_qc prof_* psal psal_qc sqc dir_name indnum temp temp_qc tqc
pres = pres2; psal = psal2; temp = temp2; clear pres2 psal2 temp2 pres_adj pres_qc_adj psal_adj psal_qc_adj temp_adj temp_qc_adj time_adj

save('/home/alison/Matlab/Mat files/Global Mapping 2/North Pacific/prof_2009_qc.mat');

clear; load /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat
save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats_2009.mat


%% quality control

% 2010 	
clear; pack
load /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/prof_2010.mat;
load /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/prof_adj_2010.mat;
load /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/prof_err_2010.mat;

pres_qc = pres_qc';
psal_qc = psal_qc';
temp_qc = temp_qc';
pres_qc_adj = pres_qc_adj';
psal_qc_adj = psal_qc_adj';
temp_qc_adj = temp_qc_adj';

longitude(longitude < 0) = longitude(longitude < 0) + 360;

% plot locations of data
figure(1); clf
axesm('mercator','MapLatLimit',[-5 80],'MapLonLimit',[100 280])
plotm(latitude,longitude,'.');
geoshow('landareas.shp','FaceColor',[.9 .9 .9])
setm(gca, 'grid','on','meridianlabel','on','parallellabel','on','plinelocation',5,'mlinelocation',5,'plabellocation',5,'mlabellocation',5)
axis tight

figure(2); clf
plot(longitude,latitude,'.'); grid on;

% lat/lon range
ind = find(longitude <= 283 & longitude >= 119.5 & latitude <= 63 & latitude >= 0);
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

% remove sea of japan profiles
ind = find(~((longitude < 138 & latitude < 44.1 & latitude > 34.5) | (longitude < 141.1 & latitude < 45.1 & latitude > 37) | (longitude < 141.7 & latitude < 45.3 & latitude > 44))); 
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

ind = find(~(longitude > 250 & latitude > 22 & latitude < 57));
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

ind = find(~(longitude > 237 & latitude > 47.5 & latitude < 49));
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

ind = find(~(longitude < 146 & latitude > 44 & latitude < 47));
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

ind = find(~(longitude < 126.05 & latitude > 10 & latitude < 10.05));
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

ind = find(~(longitude < 125 & latitude > 11.2 & latitude < 11.3));
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

ind = find(~(longitude < 120 & latitude > 18.5 & latitude < 18.8));
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

ind = find(~(longitude < 120 & latitude > 4.9 & latitude < 4.94));
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

ind = find(~(longitude < 124.4 & latitude > 13.6 & latitude < 13.7));
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

% plot locations of data
figure(1); clf
axesm('mercator','MapLatLimit',[-5 80],'MapLonLimit',[100 280])
plotm(latitude,longitude,'.');
geoshow('landareas.shp','FaceColor',[.9 .9 .9])
setm(gca, 'grid','on','meridianlabel','on','parallellabel','on','plinelocation',10,'mlinelocation',10,'plabellocation',10,'mlabellocation',10)
axis tight

% look at maximum pressures
maxp = (max(abs(pres),[],1))';

% find & remove profiles with maxp < 800 db
indp = find(maxp < 800);
per800 = 100*length(indp)/length(maxp) %#ok<*NOPTS,*NASGU>

num_800 = [fnum(indp) time(indp)];

ind = find(maxp >= 800);
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

clear ans ind indp per800 t2 

save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_800
clear num_800 maxp

% delete all Nan profiles
clear ind; k = 1; 
for i = 1:length(fnum)
    t = temp(:,i); s = psal(:,i); p = pres(:,i);
    ss = sum(~isnan(s)); tt = sum(~isnan(t)); pp = sum(~isnan(p));
    
    if ~(tt == 0 || ss == 0 || pp == 0)
        ind(k) = i; k = k+1; %#ok<NAGROW>
    end
end
clear t s p ss tt pp k i
ind = ind';
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
datamode = datamode(ind); time_adj = time_adj(ind);
temp_adj = temp_adj(:,ind); pres_adj = pres_adj(:,ind); psal_adj = psal_adj(:,ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);

indt = find(sum(~isnan(temp),2) == 0); %#ok<MXFND>
indp = find(sum(~isnan(pres),2) == 0); %#ok<MXFND>
inds = find(sum(~isnan(psal),2) == 0); %#ok<MXFND>
ind = max([min(indt) min(indp) min(inds)]);
temp = temp(1:ind,:); psal = psal(1:ind,:); pres = pres(1:ind,:);
temp_qc = temp_qc(1:ind,:); pres_qc = pres_qc(1:ind,:); psal_qc = psal_qc(1:ind,:);
temp_adj = temp_adj(1:ind,:); psal_adj = psal_adj(1:ind,:); pres_adj = pres_adj(1:ind,:);
temp_qc_adj = temp_qc_adj(1:ind,:); pres_qc_adj = pres_qc_adj(1:ind,:); psal_qc_adj = psal_qc_adj(1:ind,:);
temp_err = temp_err(1:ind,:); psal_err = psal_err(1:ind,:); pres_err = pres_err(1:ind,:);
clear indt indp inds ind

save /home/alison/Matlab/M'at files'/Global' Mapping 2'/North' Pacific'/tmp.mat

% graylist floats
cd('/home/alison/Matlab/M files/Global Mapping 2/');
[grayfloat,~,start_date,end_date,errtype] = graylist('/home/alison/Matlab/Mat files/Global Mapping 2/North Pacific/');


%% %%%% quality control flags %%%%

% examine flags on temperature data
[mm,nn] = size(temp);
tqc = NaN(size(temp));
for j = 1:nn
    for i = 1:mm
        if ~isnan(temp(i,j))
            if (strcmp(temp_qc_adj(i,j),'1') || strcmp(temp_qc_adj(i,j),'0') || strcmp(temp_qc_adj(i,j),'2') || strcmp(temp_qc_adj(i,j),' ') || double(temp_qc_adj(i,j))==0)
                tqc(i,j) = 1;
            else
                tqc(i,j) = 0;
            end
        end
    end
end

% examine flags on salinity data
[mm,nn] = size(psal);
sqc = NaN(size(psal));
for j = 1:nn
    for i = 1:mm
        if ~isnan(psal(i,j))
            if (strcmp(psal_qc_adj(i,j),'1') || strcmp(psal_qc_adj(i,j),'0') || strcmp(psal_qc_adj(i,j),'2') || strcmp(psal_qc_adj(i,j),' ') || double(psal_qc_adj(i,j))==0)
                sqc(i,j) = 1;
            else
                sqc(i,j) = 0;
            end
        end
    end
end

% examine flags on pressure data
[mm,nn] = size(pres);
pqc = NaN(size(pres));
for j = 1:nn
    for i = 1:mm
        if ~isnan(pres(i,j))
            if (strcmp(pres_qc_adj(i,j),'1') || strcmp(pres_qc_adj(i,j),'0') || strcmp(pres_qc_adj(i,j),'2') || strcmp(pres_qc_adj(i,j),' ') || double(pres_qc_adj(i,j))==0)
                pqc(i,j) = 1;
            else
                pqc(i,j) = 0;
            end
        end
    end
end

clear mm nn i j


% display number of profiles with one or more bad qc flag
num_bad_temp = length(find(sum(tqc == 0,1) > 0))
num_bad_pres = length(find(sum(pqc == 0,1) > 0))
num_bad_psal = length(find(sum(sqc == 0,1) > 0))
clear num_bad_*

% make 2nd variable set with data with qc ~= 0,1,2,' ' set to NaN
temp2 = temp_adj; temp2(tqc == 0) = NaN;
psal2 = psal_adj; psal2(sqc == 0) = NaN;
pres2 = pres_adj; pres2(pqc == 0) = NaN;

% need t, s, and p at every point - make NaN's consistent across all data
ind = find(isnan(pres2));
temp2(ind) = NaN;
psal2(ind) = NaN;

ind = find(isnan(psal2));
temp2(ind) = NaN;
pres2(ind) = NaN;

ind = find(isnan(temp2));
pres2(ind) = NaN;
psal2(ind) = NaN;

% delete profiles that are all NaNs
clear ind indb; k = 1; l = 1;
for i = 1:length(fnum)
    t = temp2(:,i); s = psal2(:,i); p = pres2(:,i);
    ss = sum(~isnan(s)); tt = sum(~isnan(t)); pp = sum(~isnan(p));
    
    if ~(tt == 0 || ss == 0 || pp == 0)
        ind(k) = i; k = k+1; %#ok<SAGROW>
    else
        indb(l) = i; l = l+1; %#ok<SAGROW>
    end
end
clear t s p ss tt pp k i l
ind = ind'; indb = indb';

num_badQC = [fnum(indb) time(indb)];
save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_badQC -APPEND
clear indb num_badQC

temp2 = temp2(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres2 = pres2(:,ind);
psal2 = psal2(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
pqc = pqc(:,ind); pres = pres(:,ind); psal = psal(:,ind); sqc = sqc(:,ind); temp = temp(:,ind); tqc = tqc(:,ind);
pres_adj = pres_adj(:,ind); pres_err = pres_err(:,ind); pres_qc = pres_qc(:,ind);
psal_adj = psal_adj(:,ind); psal_err = psal_err(:,ind); psal_qc = psal_qc(:,ind);
temp_adj = temp_adj(:,ind); temp_err = temp_err(:,ind); temp_qc = temp_qc(:,ind);
datamode = datamode(ind); time_adj = time_adj(ind);

clear ind

% find & remove profiles with maxp < 800 db
maxp = (max(abs(pres2),[],1))';
indp = find(maxp < 800);
per800 = 100*length(indp)/length(maxp) %#ok<*NOPTS,*NASGU>

num_800_2 = [fnum(indp) time(indp)];

ind = find(maxp >= 800);
temp2 = temp2(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres2 = pres2(:,ind);
psal2 = psal2(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
pqc = pqc(:,ind); pres = pres(:,ind); psal = psal(:,ind); sqc = sqc(:,ind); temp = temp(:,ind); tqc = tqc(:,ind);
pres_adj = pres_adj(:,ind); pres_err = pres_err(:,ind); pres_qc = pres_qc(:,ind);
psal_adj = psal_adj(:,ind); psal_err = psal_err(:,ind); psal_qc = psal_qc(:,ind);
temp_adj = temp_adj(:,ind); temp_err = temp_err(:,ind); temp_qc = temp_qc(:,ind);
datamode = datamode(ind); time_adj = time_adj(ind);
clear ind

clear ans ind indp per800 t2 

save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_800_2
clear num_800_2 maxp

% % compute percentage of good values in each profile
% [~,nn] = size(pres2);
% percentNaN = NaN(nn,3);
% for i = 1:nn
%     p = pqc(:,i); t = tqc(:,i); s = sqc(:,i);
%     np = sum(p == 0); nt = sum(t == 0); ns = sum(s == 0); 
%     npp = sum(~isnan(p)); ntt = sum(~isnan(t)); nss = sum(~isnan(s));
%     percentNaN(i,:) = 100*[np./npp nt/ntt ns/nss];
% end
% clear i p np lastp mm nn npp t nt ntt s ns nss
% 
% ind = find(percentNaN(:,1) > 25 | percentNaN(:,2) > 25 | percentNaN(:,3) > 25);
% for ii = 1:length(ind)
%     fnum(ind(ii))
%     figure(1); clf
%     plot(temp(:,ind(ii)),-pres(:,ind(ii)),'m.-'); hold on;
%     plot(temp2(:,ind(ii)),-pres2(:,ind(ii)),'r.-');
%     disp(temp_qc(:,ind(ii))');
%     disp(pres_qc(:,ind(ii))');
%     disp(psal_qc(:,ind(ii))');
%     figure(2); clf
%     plot(psal(:,ind(ii)),-pres(:,ind(ii)),'c.-'); hold on;
%     plot(psal2(:,ind(ii)),-pres2(:,ind(ii)),'b.-');
%     pause
% end
% clear ii ans
% 
% num_percentNaN = [fnum(ind) time(ind)];
% save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_percentNaN -APPEND
% clear num_percentNaN
% 
% temp2(:,ind) = NaN; pres2(:,ind) = NaN; psal2(:,ind) = NaN;

% find & delete short profiles
n = sum(~isnan(temp2),1)';
% for i = 1:10
%     disp(i); ind = find(n == i); length(ind)
%     figure(1); clf
%     plot(temp2(:,ind),-pres2(:,ind),'.-');
%     pause
% end

ind = find(n <= 10 & n > 0); length(ind)
num_shortprof = [fnum(ind) time(ind)];
save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_shortprof -APPEND
clear num_shortprof

temp2(:,ind) = NaN; pres2(:,ind) = NaN; psal2(:,ind) = NaN;
clear n i ind ans

% delete profiles that are all NaNs
clear ind indb; k = 1; l = 1;
for i = 1:length(fnum)
    t = temp2(:,i); s = psal2(:,i); p = pres2(:,i);
    ss = sum(~isnan(s)); tt = sum(~isnan(t)); pp = sum(~isnan(p));
    
    if ~(tt == 0 || ss == 0 || pp == 0)
        ind(k) = i; k = k+1; %#ok<SAGROW>
    else
        indb(l) = i; l = l+1; %#ok<SAGROW>
    end
end
clear t s p ss tt pp k i l
ind = ind'; indb = indb';

temp2 = temp2(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres2 = pres2(:,ind);
psal2 = psal2(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
pqc = pqc(:,ind); pres = pres(:,ind); psal = psal(:,ind); sqc = sqc(:,ind); temp = temp(:,ind); tqc = tqc(:,ind);
pres_adj = pres_adj(:,ind); pres_err = pres_err(:,ind); pres_qc = pres_qc(:,ind);
psal_adj = psal_adj(:,ind); psal_err = psal_err(:,ind); psal_qc = psal_qc(:,ind);
temp_adj = temp_adj(:,ind); temp_err = temp_err(:,ind); temp_qc = temp_qc(:,ind);
datamode = datamode(ind); time_adj = time_adj(ind);
clear ind indb 

save /home/alison/Matlab/M'at files'/Global' Mapping 2'/North' Pacific'/tmp.mat

%% remove greylist floats
k = 1;
badfloats = 0; type_err = cell(1,1); bad_start = 0; bad_end = 0;
for i = 1:length(grayfloat)
    if sum(grayfloat(i) == fnum) > 0
    if grayfloat(i) ~= badfloats
        badfloats(k,1) = grayfloat(i); %#ok<SAGROW>
        type_err(k,1) = errtype(i);
        bad_start(k,1) = start_date(i); %#ok<SAGROW>
        bad_end(k,1) = end_date(i); %#ok<SAGROW>
        k = k+1;
    end
    end
end
clear k i

% ind = find(~strcmp(type_err,'potential problem with bin assignment'));
% badfloats = badfloats(ind);
% type_err = type_err(ind);
% bad_start = bad_start(ind);
% bad_end = bad_end(ind);

badstart = cellstr(num2str(bad_start));
badend = cellstr(num2str(bad_end));

for i = 1:length(badstart)
    bs = char(badstart(i));
    bs = datenum([sscanf(bs(1:4),'%f'),sscanf(bs(5:6),'%f'),sscanf(bs(7:8),'%f')]);
    bad_st(i,1) = bs; %#ok<SAGROW>
    be = char(badend(i));
    if ~strcmp(be,'     NaN') && ~strcmp(be,'NaN')
        be = datenum([sscanf(be(1:4),'%f'),sscanf(be(5:6),'%f'),sscanf(be(7:8),'%f')]);
        bad_end(i) = be;
    else
        bad_end(i) = NaN;
    end
end

fnumb = []; timeb = [];
for j = 1:length(badfloats)
    if isnan(bad_end(j))
        ind = find(~((fnum == badfloats(j))&(time > bad_st(j))));
        indb = find((fnum == badfloats(j)) & (time > bad_st(j)));
    else
        ind = find(~(fnum == badfloats(j) & time > bad_st(j) & time < bad_end(j)));
        indb = find(fnum == badfloats(j) & time > bad_st(j) & time < bad_end(j));
    end
fnumb = [fnumb; fnum(indb)]; %#ok<AGROW>
timeb = [timeb; time(indb)]; %#ok<AGROW>

temp2 = temp2(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres2 = pres2(:,ind);
psal2 = psal2(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
pqc = pqc(:,ind); pres = pres(:,ind); psal = psal(:,ind); sqc = sqc(:,ind); temp = temp(:,ind); tqc = tqc(:,ind);
pres_adj = pres_adj(:,ind); pres_err = pres_err(:,ind); pres_qc = pres_qc(:,ind);
psal_adj = psal_adj(:,ind); psal_err = psal_err(:,ind); psal_qc = psal_qc(:,ind);
temp_adj = temp_adj(:,ind); temp_err = temp_err(:,ind); temp_qc = temp_qc(:,ind);
datamode = datamode(ind); time_adj = time_adj(ind);

end

num_gray = [fnumb timeb];
save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_gray -APPEND

clear num_gray fnumb timeb indb ind
clear bad* end_date j ind param start_date be bs i 
clear ans errtype grayfloat ii percentNaN type_err

% graylist floats
cd('/home/alison/Matlab/M files/Global Mapping 2/');
[grayfloat,param,start_date,end_date,errtype] = graylist('/home/alison/Matlab/Mat files/Global Mapping 2/North Pacific/');

save /home/alison/Matlab/M'at files'/Global' Mapping 2'/North' Pacific'/tmp.mat

%% look for out-of-range pressures
figure(1); clf; imagesc(pres); colorbar;
figure(2); clf; imagesc(pres2); colorbar

outofrangenum = [];

ind = find(pres > 2200); length(ind)
[i,j] = find(pres2 > 2200); length(i)
tempnum = fnum(j);
presnum = 0; k = 1;
for ii = 1:length(tempnum)
if tempnum(ii) ~= presnum
    presnum(k,1) = tempnum(ii); %#ok<SAGROW>
    k = k+1;
end
end
outofrangenum = [outofrangenum; presnum];

% for ii = 1:length(i)
%     figure(1); clf
%     plot(temp2(i(ii),:),-pres2(i(ii),:),'r.-'); hold on
%     plot(psal2(i(ii),:),-pres2(i(ii),:),'b.-');
%     pres2(i(ii),1:30)
%     pause
% end

ind = find(pres < 0); length(ind)
[~,j] = find(pres2 < 0); length(j)
tempnum = fnum(j);
presnum = 0; k = 1;
for ii = 1:length(tempnum)
if tempnum(ii) ~= presnum
    presnum(k,1) = tempnum(ii); %#ok<SAGROW>
    k = k+1;
end
end
outofrangenum = [outofrangenum; presnum];

% for ii = 1:length(i)
%     figure(1); clf
%     plot(temp(i(ii),:),-pres(i(ii),:),'r.-'); hold on
%     plot(psal(i(ii),:),-pres(i(ii),:),'b.-');
%     pres2(i(ii),1:30)
%     pause
% end

% look for out-of-range temperatures
figure(1); clf; imagesc(temp); colorbar;
figure(2); clf; imagesc(temp2); colorbar

ind = find(temp > 40); length(ind)
[~,j] = find(temp2 > 40); length(j)
tempnum = fnum(j);
tmpnum = 0; k = 1;
for ii = 1:length(tempnum)
if tempnum(ii) ~= tmpnum
    tmpnum(k,1) = tempnum(ii); %#ok<SAGROW>
    k = k+1;
end
end
outofrangenum = [outofrangenum; tmpnum];

% for ii = 1:length(i)
%     figure(1); clf
%     plot(temp2(i(ii),:),-pres2(i(ii),:),'r.-'); hold on
%     plot(psal2(i(ii),:),-pres2(i(ii),:),'b.-');
%     temp2(i(ii),1:30)
%     pause
% end

ind = find(temp <= -3); length(ind)
[~,j] = find(temp2 <= -3); length(j)
tempnum = fnum(j);
tmpnum = 0; k = 1;
for ii = 1:length(tempnum)
if tempnum(ii) ~= tmpnum
    tmpnum(k,1) = tempnum(ii); %#ok<SAGROW>
    k = k+1;
end
end
outofrangenum = [outofrangenum; tmpnum];

% for ii = 1:length(i)
%     figure(1); clf
%     plot(temp2(i(ii),:),-pres2(i(ii),:),'r.-'); hold on
%     plot(psal2(i(ii),:),-pres2(i(ii),:),'b.-');
%     temp2(i(ii),1:30)
%     pause
% end

% look for out-of-range salinities
figure(1); clf; imagesc(psal); colorbar;
figure(2); clf; imagesc(psal2); colorbar

ind = find(psal > 39); length(ind)
[~,j] = find(psal2 > 39); length(j)
tempnum = fnum(j);
psalnum = 0; k = 1;
for ii = 1:length(tempnum)
if tempnum(ii) ~= psalnum
    psalnum(k,1) = tempnum(ii); %#ok<SAGROW>
    k = k+1;
end
end
outofrangenum = [outofrangenum; psalnum];

% for ii = 1:length(i)
%     figure(1); clf
% %     plot(temp2(i(ii),:),-pres2(i(ii),:),'r.-'); hold on
%     plot(psal2(:,j(ii)),-pres2(:,j(ii)),'b.-'); grid on
% %     psal2(i(ii),1:30)
%     pause
% end

ind = find(psal <= 20); length(ind)
[~,j] = find(psal2 <= 20); length(j)
tempnum = fnum(j);
psalnum = 0; k = 1;
for ii = 1:length(tempnum)
if tempnum(ii) ~= psalnum
    psalnum(k,1) = tempnum(ii); %#ok<SAGROW>
    k = k+1;
end
end
outofrangenum = [outofrangenum; psalnum];

% for ii = 1:length(i)
%     figure(1); clf
% %     plot(temp2(i(ii),:),-pres2(i(ii),:),'r.-'); hold on
%     plot(psal2(:,j(ii)),-pres2(:,j(ii)),'b.-'); grid on
% %     psal2(i(ii),1:30)
%     pause
% end

badnum = sort(outofrangenum);

clear psalnum tmpnum tempnum presnum i ii ind j k outofrangenum NaNnum ans

badnum = badnum(5:end)

fnumb = []; timeb = [];

badnum(i)
ind = find(fnum == badnum(i));
figure(1); clf; plot(temp2(:,ind),-pres2(:,ind),'.-'); grid on;
figure(2); clf; plot(psal2(:,ind),-pres2(:,ind),'.-'); grid on;
tmpt = temp2(:,ind); tmps = psal2(:,ind); tmpp = pres2(:,ind); tmpj = datestr(time(ind));
if sum(badnum(i) == grayfloat) > 0
  disp('greylisted');
  ii = find(badnum(i) == grayfloat);
  disp(param(ii)); disp(errtype(ii)); disp(start_date(ii)); disp(end_date(ii));
end

i = 1;
% set bad values to Nan
[ii,j] = find(psal2(:,ind) < 22);
for k = 1:length(ii)
    pres2(ii(k),ind(j(k))) = NaN;
end
clear ii j k

i = 2;
% set bad values to Nan
temp2(1,ind(22)) = NaN; 
psal2(isnan(temp2)) = NaN;
fnumb = [fnumb; fnum(ind(22))];
timeb = [timeb; time(ind(22))];

num_lowS = [fnumb timeb];
save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_lowS -APPEND

clear num_lowS fnumb timeb indb ind ans i ii tmpj tmpp tmps tmpt

% for ii = 1:length(ind)
%     ii
%     figure(1); clf; plot(psal2(ind,:)',-pres2(ind,:)','-','color',[.9 .9 .9]);
%     hold on; plot(psal2(ind(ii),:)',-pres2(ind(ii),:)','r.-'); grid on;
%     figure(2); clf;
%     axesm('mercator','MapLatLimit',[min(latitude(ind))-5 max(latitude(ind))+5],'MapLonLimit',[min(longitude(ind))-5 max(longitude(ind))+5])
%     plotm(latitude(ind),longitude(ind),'.-'); hold on;
%     plotm(latitude(ind(ii)),longitude(ind(ii)),'r.');
%     geoshow('landareas.shp','FaceColor',[.9 .9 .9])
%     setm(gca, 'grid','on','meridianlabel','on','parallellabel','on','plinelocation',10,'mlinelocation',10,'plabellocation',10,'mlabellocation',10)
%     axis tight;
%     figure(3); clf; plot(temp2(ind,:)',-pres2(ind,:)','-','color',[.9 .9 .9]);
%     hold on; plot(temp2(ind(ii),:)',-pres2(ind(ii),:)','r.-'); grid on;
%     pause
% end

% delete profiles that are all NaNs
clear ind indb; k = 1; l = 1;
for i = 1:length(fnum)
    t = temp2(:,i); s = psal2(:,i); p = pres2(:,i);
    ss = sum(~isnan(s)); tt = sum(~isnan(t)); pp = sum(~isnan(p));
    
    if ~(tt == 0 || ss == 0 || pp == 0)
        ind(k) = i; k = k+1; %#ok<SAGROW>
    else
        indb(l) = i; l = l+1; %#ok<SAGROW>
    end
end
clear t s p ss tt pp k i l
ind = ind'; indb = indb';

temp2 = temp2(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres2 = pres2(:,ind);
psal2 = psal2(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
pqc = pqc(:,ind); pres = pres(:,ind); psal = psal(:,ind); sqc = sqc(:,ind); temp = temp(:,ind); tqc = tqc(:,ind);
pres_adj = pres_adj(:,ind); pres_err = pres_err(:,ind); pres_qc = pres_qc(:,ind);
psal_adj = psal_adj(:,ind); psal_err = psal_err(:,ind); psal_qc = psal_qc(:,ind);
temp_adj = temp_adj(:,ind); temp_err = temp_err(:,ind); temp_qc = temp_qc(:,ind);
datamode = datamode(ind); time_adj = time_adj(ind);
clear ind indb 

%% find profiles with repeat pressures
clear ind1; k = 1;
for i = 1:length(fnum)
    p = pres2(:,i);
    for j = 1:length(p)
        if ~isnan(p(j))
            if sum(p == p(j)) > 1
                ind1(k) = i; %#ok<SAGROW>
                k = k+1;
            end
        end
    end
end
clear p i j k
ind1 = ind1';

% individual float numbers of profiles with repeat pressures
k = 1; ind = 0;
for i = 1:length(ind1)
    if fnum(ind1(i)) ~= ind
        ind(k) = fnum(ind1(i)); %#ok<SAGROW>
        k = k+1;
    end
end
ind = ind'; clear i k
length(ind)

% for i = 1:length(ind)
%     ii = find(fnum == ind(i));
%     plot(temp2(:,ii),-pres2(:,ii),'.-');
%     ind(i)
%     pres2(1:2,ii)
%     pause
% end
% clear ii i

i = 1;

i = i+1;
ind(i)
ii = find(fnum == ind(i));
figure(1); clf; plot(temp2(:,ii),-pres2(:,ii),'.-'); grid on;
figure(2); clf; plot(psal2(:,ii),-pres2(:,ii),'.-'); grid on;
tmpt = temp2(:,ii); tmps = psal2(:,ii); tmpp = pres2(:,ii); tmpj = datestr(time(ii));
if sum(ind(i) == grayfloat) > 0
  disp('greylisted');
  jj = find(ind(i) == grayfloat);
  disp(param(jj)); disp(errtype(jj)); disp(start_date(jj)); disp(end_date(jj));
end

clear ind2 ind3;
k = 1;
for iii = 1:length(ii)
p = pres2(:,ii(iii));
    for j = 1:length(p)
        if ~isnan(p(j))
            if sum(p == p(j)) > 1
                ind2(k) = iii;  %#ok<SAGROW>
                ind3(k) = j; %#ok<SAGROW>
                k = k+1;
            end
        end
    end
end
disp(ind2)

temp2(ind3(1),ii(ind2(1))) = NaN;

ind = find(isnan(pres2));
temp2(ind) = NaN;
psal2(ind) = NaN;

ind = find(isnan(psal2));
temp2(ind) = NaN;
pres2(ind) = NaN;

ind = find(isnan(temp2));
pres2(ind) = NaN;
psal2(ind) = NaN;

% find profiles with increasing pressures 
k = 1; ind1 = [];
for i = 1:length(fnum)
    p = pres2(:,i);
    p = p(~isnan(p));
    pd = p(2:end)-p(1:end-1);
    if sum(pd < 0) > 0
        ind1(k) = i; %#ok<SAGROW>
        k = k+1;
    end
end
clear ind i k p pd s t j
ind1 = ind1';

% individual float numbers of profiles with increasing pressures
k = 1; ind = 0;
for i = 1:length(ind1)
    if fnum(ind1(i)) ~= ind
        ind(k) = fnum(ind1(i)); %#ok<SAGROW>
        k = k+1;
    end
end
ind = ind'; clear i k
length(ind)

for i = 1:length(ind)
    ii = find(fnum == ind(i));
    plot(temp2(:,ii),-pres2(:,ii),'.-');
    ind(i)
    pres2(1:2,ii)
    pause
end
clear ii i

i = 1;
pres2(35,ii(3)) = NaN; 

ind(i)
ii = find(fnum == ind(i));
figure(1); clf; plot(temp2(:,ii),-pres2(:,ii),'.-'); grid on;
figure(2); clf; plot(psal2(:,ii),-pres2(:,ii),'.-'); grid on;
tmpt = temp2(:,ii); tmps = psal2(:,ii); tmpp = pres2(:,ii); tmpj = datestr(time(ii));
if sum(ind(i) == grayfloat) > 0
  disp('greylisted');
  jj = find(ind(i) == grayfloat);
  disp(param(jj)); disp(errtype(jj)); disp(start_date(jj)); disp(end_date(jj));
end

k = 1; clear ind2
for iii = 1:length(ii)
    p = pres2(:,ii(iii));
    p = p(~isnan(p));
    pd = p(2:end)-p(1:end-1);
    if sum(pd < 0) > 0
        ind2(k) = iii; %#ok<SAGROW>
        k = k+1;
    end
end
disp(ind2)

ind = find(isnan(pres2));
temp2(ind) = NaN;
psal2(ind) = NaN;

ind = find(isnan(psal2));
temp2(ind) = NaN;
pres2(ind) = NaN;

ind = find(isnan(temp2));
pres2(ind) = NaN;
psal2(ind) = NaN;

clear ans i ii ind ind1 jj badnum tmpp tmpt tmpj tmps
clear iii ind2 ind3

%% plot all profiles and look for errant ones
figure(1); clf
plot(temp2,-pres2,'.-'); grid on; axis tight;

figure(2); clf 
plot(psal2,-pres2,'.-'); grid on; axis tight;

pres2(:,j) = NaN;

[i,j] = find(pres2 > 1040 & psal2 < 34.15);
for k = 1:length(i)
    pres2(i(k),j(k)) = NaN;
end

[i,j] = find(pres2 > 695 & pres2 < 700 & psal2 < 33.42);
for k = 1:length(i)
    pres2(i(k),j(k)) = NaN;
end

[i,j] = find(pres2 > 970 & pres2 < 980 & psal2 < 33.9);
for k = 1:length(i)
    pres2(i(k),j(k)) = NaN;
end

[i,j] = find(pres2 > 1350 & pres2 < 1450 & psal2 < 34.3);
for k = 1:length(i)
    pres2(i(k),j(k)) = NaN;
end

% i = 6;
% temp2(:,ind(12)) = NaN;
% psal2(isnan(temp2)) = NaN;
% fnumb = [fnumb; fnum(ind(12))]; timeb = [timeb; time(ind(12))];
% 
% badnum(i)
% ind = find(fnum == badnum(i));
% figure(1); clf; plot(temp2(:,ind),-pres2(:,ind),'.-'); grid on;
% figure(2); clf; plot(psal2(:,ind),-pres2(:,ind),'.-'); grid on;
% tmpt = temp2(:,ind); tmps = psal2(:,ind); tmpp = pres2(:,ind); tmpj = datestr(time(ind));
% if sum(badnum(i) == grayfloat) > 0
%   disp('greylisted');
%   ii = find(badnum(i) == grayfloat);
%   disp(param(ii)); disp(errtype(ii)); disp(start_date(ii)); disp(end_date(ii));
% end
% 
% num_errant = [fnumb timeb];
% save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_errant -APPEND
% clear num_errant fnumb timeb i ii ind indnum j jj badnum ans tmpj tmpp tmps tmpt


% need t, s, and p at every point - make NaN's consistent across all data
ind = find(isnan(pres2));
temp2(ind) = NaN;
psal2(ind) = NaN;

ind = find(isnan(psal2));
temp2(ind) = NaN;
pres2(ind) = NaN;

ind = find(isnan(temp2));
pres2(ind) = NaN;
psal2(ind) = NaN;

clear ind

% delete profiles that are all NaNs
clear ind indb; k = 1; l = 1;
for i = 1:length(fnum)
    t = temp2(:,i); s = psal2(:,i); p = pres2(:,i);
    ss = sum(~isnan(s)); tt = sum(~isnan(t)); pp = sum(~isnan(p));
    
    if ~(tt == 0 || ss == 0 || pp == 0)
        ind(k) = i; k = k+1; %#ok<SAGROW>
    else
        indb(l) = i; l = l+1; %#ok<SAGROW>
    end
end
clear t s p ss tt pp k i l
ind = ind'; indb = indb';

temp2 = temp2(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres2 = pres2(:,ind);
psal2 = psal2(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
pqc = pqc(:,ind); pres = pres(:,ind); psal = psal(:,ind); sqc = sqc(:,ind); temp = temp(:,ind); tqc = tqc(:,ind);
pres_adj = pres_adj(:,ind); pres_err = pres_err(:,ind); pres_qc = pres_qc(:,ind);
psal_adj = psal_adj(:,ind); psal_err = psal_err(:,ind); psal_qc = psal_qc(:,ind);
temp_adj = temp_adj(:,ind); temp_err = temp_err(:,ind); temp_qc = temp_qc(:,ind);
datamode = datamode(ind); time_adj = time_adj(ind);
clear ind indb

clear ans badnum end_date errtype grayfloat ii ind ind1 j jj kk num_* param
clear tmpj tmpp tmps tmpt start_date
% 
%% maximum pressure
maxp = (max(abs(pres2),[],1))';

% find & remove profiles with maxp < 800 db
indp = find(maxp < 800);
per800 = 100*length(indp)/length(maxp)

num_800_3 = [fnum(indp) time(indp)];

ind = find(maxp >= 800);
temp2 = temp2(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres2 = pres2(:,ind);
psal2 = psal2(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
pqc = pqc(:,ind); pres = pres(:,ind); psal = psal(:,ind); sqc = sqc(:,ind); temp = temp(:,ind); tqc = tqc(:,ind);
pres_adj = pres_adj(:,ind); pres_err = pres_err(:,ind); pres_qc = pres_qc(:,ind);
psal_adj = psal_adj(:,ind); psal_err = psal_err(:,ind); psal_qc = psal_qc(:,ind);
temp_adj = temp_adj(:,ind); temp_err = temp_err(:,ind); temp_qc = temp_qc(:,ind);
datamode = datamode(ind); time_adj = time_adj(ind);

clear ans ind indp per800 t2 

% % minimum pressure
% minp = (min(abs(pres2),[],1))';
% 
% % find & remove profiles with maxp < 800 db
% indp = find(minp > 400);
% per400 = 100*length(indp)/length(minp)
% 
% num_400 = [fnum(indp) time(indp)];
% 
% ind = find(minp <= 400);
% temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
% psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
% temp_qc = temp_qc(:,ind); pres_qc = pres_qc(:,ind); psal_qc = psal_qc(:,ind);
% position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
% prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
% temp2 = temp2(:,ind); psal2 = psal2(:,ind); pres2 = pres2(:,ind);
% tqc = tqc(:,ind); pqc = pqc(:,ind); sqc = sqc(:,ind);
% 
% clear ans ind indp per800 t2 minp maxp per400

save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat num_800_3 -APPEND
clear num_800_3 num_400 maxp

clear ind3 jjj outofrangenum psalnum presnum tempnum tmpnum

%% remove fsi sensors
ind = find((type == 842 | (type == 847 | (type == 852 | type == 857))));

k = 1; indnum = 0;
for i = 1:length(fnum(ind))
    if fnum(ind(i)) ~= indnum
        indnum(k) = fnum(ind(i)); %#ok<SAGROW>
        k = k+1;
    end
end
indnum = indnum';
clear i k

cd /home/alison/Matlab/M' files'/Global' Mapping 2'/

for i = 1:length(indnum)
    fn = num2str(indnum(i));
    
    filename = strcat('/home/alison/Data/Argo Profiles/Float/North Pacific/','Profiles/',fn,'_prof.nc');
    nc = netcdf.open(char(filename),'NC_NOWRITE');
    pi = ncchar(nc,'PI_NAME');
    PI(i,:) = (pi(:,1))';
    netcdf.close(nc);
end

indnum2 = NaN(length(indnum),1);
for i = 1:length(indnum)
   
    if strcmp(PI(i,1:5),'BRECK')
        indnum2(i) = 1;
    else indnum2(i) = 0;
    end
end

indnum = indnum(indnum2==0);

for i = 1:length(indnum)

    ind = find(fnum == indnum(i));
    time(ind) = NaN;
end

ind = find(~isnan(time));

temp2 = temp2(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres2 = pres2(:,ind);
psal2 = psal2(:,ind); fnum = fnum(ind); type = type(ind);
temp_qc_adj = temp_qc_adj(:,ind); pres_qc_adj = pres_qc_adj(:,ind); psal_qc_adj = psal_qc_adj(:,ind);
position_qc = position_qc(ind); juld_qc = juld_qc(ind); prof_temp_qc = prof_temp_qc(ind);
prof_psal_qc = prof_psal_qc(ind); prof_pres_qc = prof_pres_qc(ind);
pqc = pqc(:,ind); pres = pres(:,ind); psal = psal(:,ind); sqc = sqc(:,ind); temp = temp(:,ind); tqc = tqc(:,ind);
pres_adj = pres_adj(:,ind); pres_err = pres_err(:,ind); pres_qc = pres_qc(:,ind);
psal_adj = psal_adj(:,ind); psal_err = psal_err(:,ind); psal_qc = psal_qc(:,ind);
temp_adj = temp_adj(:,ind); temp_err = temp_err(:,ind); temp_qc = temp_qc(:,ind);
datamode = datamode(ind); time_adj = time_adj(ind);

clear PI ans filename fn i ind indnum indnum2 nc pi j


%% %%%%%%%%%%%%%%%
save('/home/alison/Matlab/Mat files/Global Mapping 2/North Pacific/prof_2010_qc_all.mat');   

clear juld_qc percentNaN position_qc pqc pres pres_qc prof_* psal psal_qc sqc dir_name indnum temp temp_qc tqc
pres = pres2; psal = psal2; temp = temp2; clear pres2 psal2 temp2 pres_adj pres_qc_adj psal_adj psal_qc_adj temp_adj temp_qc_adj time_adj

save('/home/alison/Matlab/Mat files/Global Mapping 2/North Pacific/prof_2010_qc.mat');

clear; load /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats.mat
save /home/alison/Matlab/Ma't files'/Global' Mapping 2'/North' Pacific'/bad_floats_2010.mat




%% find individual float numbers
%%%%%%%%%%%%%%%%%%
k = 1; indnum = 0;
for i = 1:length(fnum)
    if fnum(i) ~= indnum
        indnum(k) = fnum(i); %#ok<SAGROW>
        k = k+1;
    end
end
indnum = indnum';
%%%%%%%%%%%%%%%%%%
clear i k
