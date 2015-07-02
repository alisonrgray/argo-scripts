% 
% Save profile adjusted data to .mat files (do a few thousand profiles at a time for memory reasons)
%
% Inputs are name of the directory holding all the argo profile data
% (example: /Users/henridrake/princeton/Data/Argo/Argo_profs/)
%
function [temp_adj,pres_adj,psal_adj,time_adj,datamode]=save_profiles_adj(dir_name,num_files)

% For this function:
% 1) From profiles, need to create variables.
% 2) Clear all variables except the ones we want to save.

% Go through the directory and pick out the profiles we want.
% strtok takes off the _prof.nc part. We only want the WMO number.
% Also, take out the first 2 files in the directory which are . and ..

cd(dir_name);
direc = dir;
floats = {};
[floats{1:length(direc),1}] = deal(direc.name);
floats = strtok(floats,'_');
floats = floats(3:end);

total_floats=zeros(length(floats));
for i = 1:length(floats);
    total_floats(i) = str2double(char(floats(i)));
end
clear direc floats i
total_numfloats = length(total_floats);
disp(total_numfloats)
step=ceil(total_numfloats/num_files);
disp((num_files-1)*step+1)

% Iterate through all files... sum of all the sizes of the _prof.nc files.
for sp = 1:num_files
    if sp==num_files
        some_floats=total_floats((num_files-1)*step+1:total_numfloats);
        numfloats=length(some_floats);
    else
        some_floats=total_floats((sp-1)*step+1:sp*step);
        numfloats=length(some_floats);
    end

    i=1;
    
    % This for loop adds up all the sizes of variables p for each profile.
    % What is p? The number of profiles per float maybe?
    for n = 1:numfloats
        
        fn = num2str(some_floats(n));
        
        % disp(fn)
        if fn=='0'
            break
        end
        
        filename = strcat(dir_name,fn,'_prof.nc');

        % Read each profile
        cd('/Users/henridrake/Scripts/Matlab/argo_scripts/')
        [~,p] = read_prof(filename);

        % What is this for? Still don't know
        if some_floats(n) == 2900171
            p = [p(:,1:28) p(:,30:end)]; p = p(1:50,:);
        end


        [~,nn] = size(p);
        %disp(qq)
        % What is p?
       
        i = i+nn;

        clear p nn mm filename

    end

    % Max number of measurements per profile = 7200
    temp_adj = NaN(7200,i); pres_adj = NaN(7200,i); psal_adj = NaN(7200,i);
    datamode = char(NaN(i,1)); time_adj = NaN(i,1);
    
    i = 1;
    for n = 1:numfloats

        fn = num2str(some_floats(n));
        
        if fn=='0'
            break
        end
        filename = strcat(dir_name,fn,'_prof.nc');
        
        cd('/Users/henridrake/Scripts/Matlab/argo_scripts/')
        [t,p,s,~,~,~,jd,dm] = read_prof_adj(filename);
    
        % Save the size of p as [mm,nn]. What is the first column
        % if second is number of profiles?
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
        time_adj(i:i-1+nn) = jd;

        i = i+nn;

        clear t p s jd mm nn ii jj ind dm

    end

    savename=strcat('/Users/henridrake/Scripts/Matlab/mat/',num2str(sp),'_adj.mat')
    save(savename)
end

return