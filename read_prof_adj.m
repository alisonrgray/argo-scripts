% [temp,pres,psal,lon,lat,juld,type] = read_prof(filename)
% Read variables from one netcdf profile file

function [temp,pres,psal,temperr,preserr,psalerr,juld,datamode] = read_prof_adj(filename)

% open netcdf file & get number of variables
nc = netcdf.open(char(filename),'NC_NOWRITE');
[~,nvars] = netcdf.inq(nc);

% check to see if psal, temp, and pres are variables of this file
sflag = 0; tflag = 0; pflag = 0; 
for i = 1:nvars
    if strcmp(netcdf.inqVar(nc,i-1),'PSAL_ADJUSTED')
    sflag = 1;
    end
    if strcmp(netcdf.inqVar(nc,i-1),'TEMP_ADJUSTED')
    tflag = 1;
    end
    if strcmp(netcdf.inqVar(nc,i-1),'PRES_ADJUSTED')
    pflag = 1;
    end
end

% load variables
if sflag
    psal = ncread(filename,'PSAL_ADJUSTED');
    psalerr = ncread(filename,'PSAL_ADJUSTED_ERROR');
else
    psal = [];
    psalerr = [];
end

if tflag
    temp = ncread(filename,'TEMP_ADJUSTED');
    temperr = ncread(filename,'TEMP_ADJUSTED_ERROR');
else
    temp = [];
    temperr = [];
end

if pflag
    pres = ncread(filename,'PRES_ADJUSTED');
    preserr = ncread(filename,'PRES_ADJUSTED_ERROR');
else
    pres = [];
    preserr = [];
end

datamode = ncread(filename,'DATA_MODE');

juld = ncread(filename,'JULD');    

% reference date
varid = netcdf.inqVarID(nc,'REFERENCE_DATE_TIME');
reference_date_time = char(netcdf.getVar(nc,varid))';

if ~isempty(reference_date_time)
dayref=datenum(sscanf(reference_date_time(1:4),'%f'),sscanf(reference_date_time(5:6),'%f'),sscanf(reference_date_time(7:8),'%f'),0,0,0);
juld = juld+dayref;
else
    disp('no reference date, using 1950-01-01')
    juld = juld+712224;
end

netcdf.close(nc);

clear nc reference_date_time dayref varid
return
