% [temp,pres,psal,lon,lat,juld,type] = read_prof(filename)
% Read variables from one netcdf profile file

function [tqc,pqc,sqc,posqc,juldqc,proftqc,profsqc,profpqc,tqc_adj,pqc_adj,sqc_adj] = read_prof_qc(filename)

% open netcdf file & get number of variables
nc = netcdf.open(char(filename),'NC_NOWRITE');
[ndims,nvars] = netcdf.inq(nc);

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
    varid = netcdf.inqVarID(nc,'PSAL_QC');
    sqc = netcdf.getVar(nc,varid);
    varid = netcdf.inqVarID(nc,'PROFILE_PSAL_QC');
    profsqc = netcdf.getVar(nc,varid);
    varid = netcdf.inqVarID(nc,'PSAL_ADJUSTED_QC');
    sqc_adj = netcdf.getVar(nc,varid);
else
    sqc = [];
    profsqc = [];
    sqc_adj = [];
end

if tflag
    varid = netcdf.inqVarID(nc,'TEMP_QC');
    tqc = netcdf.getVar(nc,varid);
    varid = netcdf.inqVarID(nc,'PROFILE_TEMP_QC');
    proftqc = netcdf.getVar(nc,varid);
    varid = netcdf.inqVarID(nc,'TEMP_ADJUSTED_QC');
    tqc_adj = netcdf.getVar(nc,varid);
else
    tqc = [];
    proftqc = [];
    tqc_adj = [];
end

if pflag
    varid = netcdf.inqVarID(nc,'PRES_QC');
    pqc = netcdf.getVar(nc,varid);
    varid = netcdf.inqVarID(nc,'PROFILE_PRES_QC');
    profpqc = netcdf.getVar(nc,varid);
    varid = netcdf.inqVarID(nc,'PRES_ADJUSTED_QC');
    pqc_adj = netcdf.getVar(nc,varid);
else
    pqc = [];
    profpqc = [];
    pqc_adj = [];
end

varid = netcdf.inqVarID(nc,'JULD_QC');
juldqc = netcdf.getVar(nc,varid);

varid = netcdf.inqVarID(nc,'POSITION_QC');
posqc = netcdf.getVar(nc,varid);

netcdf.close(nc);

clear nc varid
return
