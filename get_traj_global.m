% Function to download list of traj.nc and associated meta.nc files
%
% 
% # Title : Trajectory directory file of the Argo global Data Assembly Center
% # Description : The directory file describes all trajectory files of the argo GDAC ftp site.
% # Project : ARGO
% # Format version : 2.0
% # Date of update : 20070924230635
% # FTP root number 1 : ftp://ftp.ifremer.fr/ifremer/argo/dac
% # FTP root number 2 : ftp://usgodae.usgodae.org/pub/outgoing/argo/dac
% # GDAC node : CORIOLIS
% # file,latitude_max,latitude_min,longitude_max,longitude_min,profiler_type,institution,date_update
%
% INPUTS: 
% target = directory name on local machine to download files to; 
% update = download all files that have been updated since this date, in format yyyymmdd000000 (enter 0 for all files);
% email = email address (for IFREMER login);
% getlist_flag = set to 1 if need to download index file; 
% fsi_flag = set to = 1 if want to disregard all provor and fsi floats
%
% OUTPUTS: 
% float_num = list of all float numbers (WMO IDs);
% float_type = list of floats that had no prof file in directory
%%%%%%%%%%%%%%%%%%
% 
% See get_prof_global.m for a better commented version of this script (they
% are almost identical)


function [float_num,float_type] = get_traj_global(target,update,email,getlist_flag,fsi_flag)

traj_file = 'ar_index_global_traj.txt';

% ftp to coriolis site and download traj_file
ff = ftp('ftp.ifremer.fr','anonymous','hdrake@princeton.edu');
cd(ff,'ifremer/argo')

% For some reason, need to be in passive mode (google matlab ftp passive if
% unsure)
pasv(ff);
 
if getlist_flag == 1
    mget(ff,traj_file);
    
    % The first non-commented line will probably be text, explaining what
    % each column is. Need to comment this out. Also, search text for other
    % signs of text. I found 3 or 4 profiles has 'null' for their type.
    % Comment those lines out.
    disp('Please check ar_index_global_prof.txt for NaN''s and text.  Put ''#'' in front of unnecessary lines.');
    disp('Please delete floats with float_types that are ''null'' in ar_index_global_traj.txt.') 
    disp('Press any key to resume');
    pause
end

fid = fopen(traj_file);
C = textscan(fid,'%s %f %f %f %f %f %s %f','Delimiter',',','CommentStyle','#');

filelist = C{1,1};
float_type = C{1,6};
date_update = C{1,8};

fclose(fid);
clear C fid

% Remove provor, fsi sensor floats
if fsi_flag == 1
    ind = find(~(float_type == 840 | (float_type == 841 | (float_type == 842 | (float_type == 847 | (float_type == 852 | float_type == 857))))));
    filelist = filelist(ind);
    date_update = date_update(ind);
end

% find floats that need to be updated
ind = find(date_update >= update);
filelist = filelist(ind);
clear date_update

disp(length(filelist))

% directory to download to.
cd(target)
float_num = NaN(size(ind));

cd(ff,'dac')

for i = 1:length(ind);

   C = textscan(char(filelist(i)),'%s %f %s','Delimiter','/');
   float_num(i,1) = double(C{1,2});
   cd(ff,strcat(char(C{1,1}),'/',num2str(C{1,2})));
   cd(strcat(target,'Trajectories/'));
   mget(ff,char(C{1,3}))
   cd(strcat(target,'Metadata/'));
   mget(ff,strcat(num2str(C{1,2}),'_meta.nc'))
   cd(ff,'..'); cd(ff,'..');
   
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
