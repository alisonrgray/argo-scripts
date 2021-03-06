%%%%%%%%%%%%%%%%%%
% Function to download list of prof.nc and associated files
%
% INPUTS: 
% target = directory name on local machine to download files to; 
% ln = longitude range [west east];
% lt = latitude range [south north];
% update = download all files that have been updated since this date, in format yyyymmdd000000 (enter 0 for all files);
% email = email address (for IFREMER login);
% getlist_flag = set to 1 if need to download index file; 
% fsi_flag = set to = 1 if want to disregard all provor and fsi floats
%
% OUTPUTS: 
% float_num = list of all float numbers (WMO IDs);
% no_prof = list of floats that had no prof file in directory
%%%%%%%%%%%%%%%%%%


function [float_num,no_prof] = get_prof(target,ln,lt,update,email,getlist_flag,fsi_flag)


prof_file = 'ar_index_global_prof.txt';

% ftp to coriolis site and download prof_file
ff = ftp('ftp.ifremer.fr','anonymous',email);
cd(ff,'ifremer/argo')
pasv(ff);

if getlist_flag == 1
    mget(ff,prof_file);
    
    disp('Please check ar_index_global_prof.txt for NaN''s and text.  Put ''#'' in front of unnecessary lines.');
    disp('Press any key to resume');
    pause
end

% read variables out of prof_file.txt
fid = fopen(prof_file);
C = textscan(fid,'%s %f %f %f %s %f %s %f','Delimiter',',','CommentStyle','#');

% String of file names
filelist = C{1,1};

% Floating points for lats and lons
lat = C{1,3};
lon = C{1,4};

% What is C{1,5}?

% String for type of float
float_type = C{1,6};

% What is C{1,7}?

% Float for date
date_update = C{1,8};

fclose(fid);
clear C fid

% convert longtitude to 360 scale
if ln(1) >= 180 || ln(2) >= 180
ind = find(lon < 0);
lon(ind) = lon(ind) + 360;
end

% remove provor, fsi sensor floats
if fsi_flag == 1
    ind = find(~(float_type == 840 | (float_type == 841 | (float_type == 842 | (float_type == 847 | (float_type == 852 | float_type == 857))))));
    filelist = filelist(ind);
    lat = lat(ind);
    lon = lon(ind);
    date_update = date_update(ind);
end
clear float_type

% Probably don't need this
% find floats in given latitude/longitude range
ind = find(lat <= lt(2) & lat >= lt(1) & (lon <= ln(2) & lon >= ln(1)));
filelist = filelist(ind);
date_update = date_update(ind);
clear lat lon

% find floats that need to be updated
filelist = filelist(date_update >= update);
clear date_update ind

% directory to download to.
cd(target)

% make a list of distinct floats
j = 1;
float_num = NaN(size(filelist));
for i = 1:length(filelist)
    C = textscan(char(filelist(i)),'%s %f %s %s','Delimiter','/');
    if i == 1
        float_num(1) = C{1,2};
        dac(1,1) = C{1,1};
        j = j+1;
    else if float_num(1:j-1) ~= C{1,2}
            float_num(j) = C{1,2};
            dac(j,1) = C{1,1}; %#ok<AGROW>
            j = j+1;
        end
    end
end
float_num = float_num(~isnan(float_num));

% display number of floats
disp(length(float_num))

% download _prof.nc file for each float
cd(ff,'dac')
k = 1;

for i = 1:length(float_num);
    cd(ff,strcat(char(dac(i)),'/',num2str(float_num(i))));
    d = dir(ff,strcat('/ifremer/argo/dac/',char(dac(i)),'/',num2str(float_num(i))));
    D = struct2cell(d); D = D(1,:);
    if sum(strcmp(D,strcat(num2str(float_num(i)),'_prof.nc')))
        mget(ff,strcat(num2str(float_num(i)),'_prof.nc'))
    else
        no_prof(k) = float_num(i); %#ok<AGROW>
        k = k+1;
    end
    cd(ff,'..'); cd(ff,'..');
    
    % this prevents issues with the connection timing out by reconnecting after every 5th file 
    if rem(i,5)==0
       close(ff);
       ff = ftp('ftp.ifremer.fr','anonymous',email);
       cd(ff,'ifremer/argo')
       pasv(ff);
       cd(ff,'dac')
    end
end

close(ff)

return
