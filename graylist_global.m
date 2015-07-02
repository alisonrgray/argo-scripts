% dirname = '/home/alison/Research/Data/Argo Profiles/Float/North Pacific/Trajectories/'
% dirname = '/home/alison/Research/Data/Argo Profiles/Float/North Pacific/Profiles/'

% flag = 1;         % to get latest greylist file from argo server
% flag = 0;         % to skip update of greylist file 

function [badnum,param,start_date,end_date,errtype] = graylist_global

grey_file = 'ar_greylist.txt';
cd('/home/alison/Data/Argo Profiles/')
% 
% if flag 
% % ftp to coriolis site and download prof_file
% ff = ftp('ftp.ifremer.fr','anonymous','alison@ocean.washington.edu');
% cd(ff,'ifremer/argo')
% pasv(ff);
% mget(ff,grey_file);
%  
% disp('Please check ar_index_greylist.txt for NaN''s and text.  Put ''#'' in front of unnecessary lines.');
% disp('Press any key to resume');
% pause
% end

fid = fopen(grey_file);
C = textscan(fid,'%f %s %f %f %f %s %s','Delimiter',',','CommentStyle','#');

badnum = C{1,1};
param = C{1,2};
start_date = C{1,3};
end_date = C{1,4};
errtype = C{1,6};

fclose(fid);
clear C fid

return



