% 
% Save profile data to .mat files
% This function splits the profiles into num_files many .mat files
%
% dir_name is the name of the directory holding all the profiles
% (example: /Users/henridrake/princeton/Data/Argo/Argo_profs/)
%
% See function get_profs_global.m to download all the profiles 
%
function [temp,pres,psal,position,juld,prof_temp,prof_psal,prof_pres]=save_profiles(dir_name,num_files)


% Go through the directory and pick out the profiles we want.
% strtok takes off the _prof.nc part. We only want the WMO number.
% Also, take out the first 2 files in the directory which are . and ..
cd(dir_name);
direc = dir;
floats = {};
[floats{1:length(direc),1}] = deal(direc.name);
floats = strtok(floats,'_');
floats = floats(3:end);

% Create an array of the float numbers
total_floats=zeros(length(floats));
for i = 1:length(floats);
    total_floats(i) = str2double(char(floats(i)));
end
clear direc floats i
total_numfloats = length(total_floats);
disp(total_numfloats)
step=ceil(total_numfloats/num_files);
disp((num_files-1)*step+1)

% Iterate through all files... sum of all the sizes of the _prof.nc files
% step is the number of profiles in each of the num_files many .mat files
for sp = 1:num_files
    
    % The last step will have less profiles because of the way we defined
    % step with the ceil() function
    % some_floats is the subset of all of the float numbers corresponding
    % to a given .mat file
    if sp==num_files
        some_floats=total_floats((num_files-1)*step+1:total_numfloats);
        numfloats=length(some_floats);
    else
        some_floats=total_floats((sp-1)*step+1:sp*step);
        numfloats=length(some_floats);
    end

    i=1;
    
    % This for loop adds up all the sizes of variables p for each profile.
    for n = 1:numfloats
        fn = num2str(some_floats(n));
        
        % For whatever reason, there are some floats missing and so the
        % zeros array does not get filled... not sure why.
        if fn=='0'
            break
        end
        
        filename = strcat(dir_name,fn,'_prof.nc');

        % Read each profile
        cd('/Users/henridrake/Scripts/Matlab/argo_scripts/')
        [~,p] = read_prof(filename);

        % What is this for? I have no idea - it was in Allison's code.
        % My guess is there is there is something weird about that float
        if some_floats(n) == 2900171
            p = [p(:,1:28) p(:,30:end)]; p = p(1:50,:);
        end

        % nn is the number of profiles for this particlar float
        [~,nn] = size(p);
       
        i = i+nn;

        clear p nn mm filename

    end

    % 7200 is the max number of measurements per profile
    % i is the total number of profiles
    % Note: the majority of floats have less than 100 measurements per
    % profile so a lot of space is wasted... not sure what to do about it
  
    longitude = NaN(i,1); latitude = NaN(i,1);  time = NaN(i,1); 
    temp = NaN(7200,i); pres = NaN(7200,i); psal = NaN(7200,i); 
    fnum = NaN(i,1); type = NaN(i,1);

    i = 1;
    for n = 1:numfloats

        fn = num2str(some_floats(n));
        
        if fn=='0'
            break
        end
        filename = strcat(dir_name,fn,'_prof.nc');
        
        cd('/Users/henridrake/Scripts/Matlab/argo_scripts/')
        [t,p,s,ln,lt,jd,ty] = read_prof(filename);

        % Why?
        if some_floats(n) == 2900171
            t = [t(:,1:28) t(:,30:end)]; t = t(1:50,:);
            s = [s(:,1:28) s(:,30:end)]; s = s(1:50,:);
            p = [p(:,1:28) p(:,30:end)]; p = p(1:50,:);
            lt = [lt(1:28); lt(30:end)];
            ln = [ln(1:28); ln(30:end)];
            jd = [jd(1:28); jd(30:end)];
            ty = [ty(1:28); ty(30:end)];
        end
    
        % Save the size of p as [mm,nn].
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

        fnum(i:i-1+nn) = some_floats(n)*ones(nn,1);

        if length(ty) == length(jd)
        type(i:i-1+nn) = ty;
        else if length(ty) == 3
                type(i:i-1+nn) = (str2num(strcat(num2str(ty(1)),num2str(ty(2)),num2str(ty(3)))))*ones(nn,1); %#ok<ST2NM>
            end
        end
        i = i+nn;

        clear t p s ln lt jd ty mm nn ii jj ind

    end

    savename=strcat('/Users/henridrake/Scripts/Matlab/mat/',num2str(sp),'.mat');
    disp(strcat('Saving ',savename))
    save(savename,'-v7.3')
end

return


