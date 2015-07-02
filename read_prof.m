% [temp,pres,psal,lon,lat,juld,type] = read_prof(filename)
% Read variables from one netcdf profile file

function [temp,pres,psal,lon,lat,juld,type] = read_prof(filename)

% open netcdf file & get number of variables
nc = netcdf.open(char(filename),'NC_NOWRITE');
% disp(netcdf.inq(nc))
[~,nvars] = netcdf.inq(nc);

% check to see if psal, temp, and pres are variables of this file
sflag = 0; tflag = 0; pflag = 0; 
for i = 1:nvars
    if strcmp(netcdf.inqVar(nc,i-1),'PSAL')
    sflag = 1;
    end
    if strcmp(netcdf.inqVar(nc,i-1),'TEMP')
    tflag = 1;
    end
    if strcmp(netcdf.inqVar(nc,i-1),'PRES')
    pflag = 1;
    end
end

% load variables
if sflag
    psal = ncread(filename,'PSAL');
else
    psal = [];
end

if tflag
    temp = ncread(filename,'TEMP');
else
    temp = [];
end

if pflag
    pres = ncread(filename,'PRES');
else
    pres = [];
end

juld = ncread(filename,'JULD');    
lat = ncread(filename,'LATITUDE'); 
lon = ncread(filename,'LONGITUDE');

% What does this do?
varid = netcdf.inqVarID(nc,'WMO_INST_TYPE');
type = str2num(char(cellstr((netcdf.getVar(nc,varid))')));  %#ok<ST2NM>

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

clear nc reference_date_time dayref varid;
return
