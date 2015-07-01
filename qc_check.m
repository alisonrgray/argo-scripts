illegal_lats_min = [];
illegal_lats_max = [];
illegal_lons_min = [];
illegal_lons_max = [];

for i = 1:50
% Load file to check out
disp(strcat('Loading: ',num2str(i),'.mat'))
load(strcat('/Users/henridrake/Scripts/Matlab/mat/',num2str(i),'.mat'))

% remove gulf of mexico
ind = find(~(longitude >= -97 & longitude <= -79 & latitude >= 19.5 & latitude <= 31));
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);

% remove sea of japan profiles
ind = find(~(longitude >= 128.7 & longitude <= 142.3 & latitude >= 35.21 & latitude <= 51.75)); 
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);

% remove mediterranean profiles
ind = find(~(longitude >= -4.208 & longitude <= 41.29 & latitude >= 30.76 & latitude <= 43.01));
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);
% remove Greenland / Norweigan Sea profiles

ind = find(~(longitude >= -23.11 & longitude <= 25.012 & latitude >= 65.65 & latitude <= 80.63));
temp = temp(:,ind); latitude = latitude(ind); longitude = longitude(ind); time = time(ind); pres = pres(:,ind);
psal = psal(:,ind); fnum = fnum(ind); type = type(ind);

% plot locations of data
figure(1); clf
axesm('mercator','MapLatLimit',[-80 80],'MapLonLimit',[-280 80])
plotm(latitude,longitude,'.');
geoshow('landareas.shp','FaceColor',[.9 .9 .9])
setm(gca, 'grid','on','meridianlabel','on','parallellabel','on','plinelocation',10,'mlinelocation',10,'plabellocation',10,'mlabellocation',10)
axis tight
[lat1,lon1]=inputm(1)
[lat2,lon2]=inputm(1)
if ~(isempty(lat1) || isempty(lat2) || isempty (lon1) || isempty (lon2))
    illegal_lats_min = horzcat(illegal_lats_min,lat1);
    illegal_lats_max = horzcat(illegal_lats_max,lat2);
    illegal_lons_min = horzcat(illegal_lons_min,lon1);
    illegal_lons_max = horzcat(illegal_lons_max,lon2);
end
   
disp('Press any key to resume');
pause
disp('Resuming')

end

save('/Users/henridrake/Scripts/Matlab/argo_scripts/marginal_seas_latslons.mat','illegal_lats_min','illegal_lats_max','illegal_lons_min','illegal_lats_max')
