% 
% Save profile data to .mat files (do a few thousand profiles at a time for memory reasons)
%
% Inputs are name of the directory holding all the argo profile data
% (example: /Users/henridrake/princeton/Data/Argo/Argo_profs/)
%
function [temp_qc,pres_qc,psal_qc,position_qc,juld_qc,prof_temp_qc,prof_psal_qc,prof_pres_qc]=save_profiles_qc(dir_name,num_files)

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
% i is our total size.
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
    position_qc = NaN(i,1); juld_qc = char(NaN(i,1)); 
    temp_qc = char(NaN(7200,i)); pres_qc = char(NaN(7200,i)); psal_qc = char(NaN(7200,i)); 
    prof_temp_qc = char(NaN(i,1)); prof_psal_qc = char(NaN(i,1));
    prof_pres_qc = char(NaN(i,1));
    temp_qc_adj = char(NaN(i,7200)); pres_qc_adj = char(NaN(i,7200)); psal_qc_adj = char(NaN(i,7200));

    i = 1;
    for n = 1:numfloats

        fn = num2str(some_floats(n));
        
        if fn=='0'
            break
        end
        
        filename = strcat(dir_name,fn,'_prof.nc');
        
        cd('/Users/henridrake/Scripts/Matlab/argo_scripts/')
        
        [tqc,pqc,sqc,posqc,juldqc,proftqc,profsqc,profpqc,tqc_adj,pqc_adj,sqc_adj] = read_prof_qc(filename);
        tqc = tqc'; pqc = pqc'; sqc = sqc';
        tqc_adj = tqc_adj'; pqc_adj = pqc_adj'; sqc_adj = sqc_adj';
        
        % Why?
        if some_floats(n) == 2900171
            tqc = [tqc(:,1:28) tqc(:,30:end)]; tqc = tqc(1:50,:);
            sqc = [sqc(:,1:28) sqc(:,30:end)]; sqc = sqc(1:50,:);
            pqc = [pqc(:,1:28) pqc(:,30:end)]; pqc = pqc(1:50,:);
            posqc = [posqc(1:28); posqc(30:end)];
            juldqc = [juldqc(1:28); juldqc(30:end)];
            proftqc = [proftqc(1:28); proftqc(30:end)];
            profsqc = [profsqc(1:28); profsqc(30:end)];
            profpqc = [profpqc(1:28); profpqc(30:end)];
        end
    
        % Save the size of p as [mm,nn]. What is the first column
        % if second is number of profiles?
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
       
        i = i+nn;

        clear tqc pqc sqc tqc_adj pqc_adj sqc_adj proftqc profpqc profsqc mm nn

    end

    savename=strcat('/Users/henridrake/Scripts/Matlab/mat/',num2str(sp),'_qc.mat');
    disp(strcat('Saving ',savename))
    save(savename,'-v7.3')
end

return

