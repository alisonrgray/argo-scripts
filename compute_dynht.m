%% compute dynamic height values

clear; pack;

cd /home/alison/Matlab/M' files'/Global' Mapping 2'/

load /home/alison/Matlab/Mat' files'/Global' Mapping 2'/North' Atlantic'/profiles_final.mat

reflevel = 900;
level = [5;10;20;30;50;75;100;125;150;200;250;300;400;500;600;700;800;900;1000;1100;1200;1300;1400;1500;1600;1700;1800;1900;2000];
lev = [0 5 10 20 30 50 75 100 125 150 175 200 250 275 300 350 400 450 500 550 600 650 700 750 800 850 900 950 1000 1050 1100 1150 1200 1250 1300 1350 1400 1450 1500 1550 1600 1650 1700 1750 1800 1850 1900 1950 2000]';

maxp = NaN(size(fnum));
minp = NaN(size(fnum));

for i = 1:length(fnum)
    
    disp(i);
    
    p = pres(:,i);
    
    mxp = max(p); mnp = min(p);
    nn = ~isnan(p);
    
    p = p(nn);
    itop = find((mnp - lev) > 0,1,'last');
    if isempty(itop); itop = 1; end
    
    ibottom = find((mxp - lev) < 0,1,'first');
    if isempty(ibottom); ibottom = length(lev); end
    
    if (mnp - lev(itop)) > .1*(mnp + lev(itop))/2  || (mnp - lev(itop)) > 4
        itop = itop+1;
    end
    
    jjj = find(p < mxp,1,'last');
    
    if (mxp - lev(ibottom)) < -.2*(mxp + lev(ibottom))/2 || (mxp - lev(ibottom)) <= -28 || lev(ibottom)-mxp >= 2*(mxp - p(jjj))
        ibottom = ibottom-1;
    end
    
    maxp(i) = lev(ibottom);
    minp(i) = lev(itop);
    
end

clear i ibottom itop p i jjj mnp mxp nn p lev

D = dynht_cont(temp,pres,psal,level,reflevel,maxp,minp);

clear indnum maxp minp pres psal temp type reflevel

save /home/alison/Matlab/Mat' files'/Global' Mapping 2'/North' Atlantic'/dynht_NA.mat

%% compute dynamic height values

clear; pack;

cd /home/alison/Matlab/M' files'/Global' Mapping 2'/

load /home/alison/Matlab/Mat' files'/Global' Mapping 2'/South' Atlantic'/profiles_final.mat

reflevel = 900;
level = [5;10;20;30;50;75;100;125;150;200;250;300;400;500;600;700;800;900;1000;1100;1200;1300;1400;1500;1600;1700;1800;1900;2000];
lev = [0 5 10 20 30 50 75 100 125 150 175 200 250 275 300 350 400 450 500 550 600 650 700 750 800 850 900 950 1000 1050 1100 1150 1200 1250 1300 1350 1400 1450 1500 1550 1600 1650 1700 1750 1800 1850 1900 1950 2000]';

maxp = NaN(size(fnum));
minp = NaN(size(fnum));

for i = 1:length(fnum)
    
    disp(i);
    
    p = pres(:,i);
    
    mxp = max(p); mnp = min(p);
    nn = ~isnan(p);
    
    p = p(nn);
    itop = find((mnp - lev) > 0,1,'last');
    if isempty(itop); itop = 1; end
    
    ibottom = find((mxp - lev) < 0,1,'first');
    if isempty(ibottom); ibottom = length(lev); end
    
    if (mnp - lev(itop)) > .1*(mnp + lev(itop))/2  || (mnp - lev(itop)) > 4
        itop = itop+1;
    end
    
    jjj = find(p < mxp,1,'last');
    
    if (mxp - lev(ibottom)) < -.2*(mxp + lev(ibottom))/2 || (mxp - lev(ibottom)) <= -28 || lev(ibottom)-mxp >= 2*(mxp - p(jjj))
        ibottom = ibottom-1;
    end
    
    maxp(i) = lev(ibottom);
    minp(i) = lev(itop);
    
end

clear i ibottom itop p i jjj mnp mxp nn p lev

D = dynht_cont(temp,pres,psal,level,reflevel,maxp,minp);

clear indnum maxp minp pres psal temp type reflevel

save /home/alison/Matlab/Mat' files'/Global' Mapping 2'/South' Atlantic'/dynht_SA.mat


%% compute dynamic height values

clear; pack;

cd /home/alison/Matlab/M' files'/Global' Mapping 2'/

load /home/alison/Matlab/Mat' files'/Global' Mapping 2'/Indian/profiles_final.mat

reflevel = 900;
level = [5;10;20;30;50;75;100;125;150;200;250;300;400;500;600;700;800;900;1000;1100;1200;1300;1400;1500;1600;1700;1800;1900;2000];
lev = [0 5 10 20 30 50 75 100 125 150 175 200 250 275 300 350 400 450 500 550 600 650 700 750 800 850 900 950 1000 1050 1100 1150 1200 1250 1300 1350 1400 1450 1500 1550 1600 1650 1700 1750 1800 1850 1900 1950 2000]';

maxp = NaN(size(fnum));
minp = NaN(size(fnum));

for i = 1:length(fnum)
    
    disp(i);
    
    p = pres(:,i);
    
    mxp = max(p); mnp = min(p);
    nn = ~isnan(p);
    
    p = p(nn);
    itop = find((mnp - lev) > 0,1,'last');
    if isempty(itop); itop = 1; end
    
    ibottom = find((mxp - lev) < 0,1,'first');
    if isempty(ibottom); ibottom = length(lev); end
    
    if (mnp - lev(itop)) > .1*(mnp + lev(itop))/2  || (mnp - lev(itop)) > 4
        itop = itop+1;
    end
    
    jjj = find(p < mxp,1,'last');
    
    if (mxp - lev(ibottom)) < -.2*(mxp + lev(ibottom))/2 || (mxp - lev(ibottom)) <= -28 || lev(ibottom)-mxp >= 2*(mxp - p(jjj))
        ibottom = ibottom-1;
    end
    
    maxp(i) = lev(ibottom);
    minp(i) = lev(itop);
    
end

clear i ibottom itop p i jjj mnp mxp nn p lev

D = dynht_cont(temp,pres,psal,level,reflevel,maxp,minp);

clear indnum maxp minp pres psal temp type reflevel

save /home/alison/Matlab/Mat' files'/Global' Mapping 2'/Indian/dynht_IND.mat




%% compute dynamic height values

clear; pack;

cd /home/alison/Matlab/M' files'/Global' Mapping 2'/

load /home/alison/Matlab/Mat' files'/Global' Mapping 2'/North' Pacific'/profiles_final.mat

reflevel = 900;
level = [5;10;20;30;50;75;100;125;150;200;250;300;400;500;600;700;800;900;1000;1100;1200;1300;1400;1500;1600;1700;1800;1900;2000];
lev = [0 5 10 20 30 50 75 100 125 150 175 200 250 275 300 350 400 450 500 550 600 650 700 750 800 850 900 950 1000 1050 1100 1150 1200 1250 1300 1350 1400 1450 1500 1550 1600 1650 1700 1750 1800 1850 1900 1950 2000]';

maxp = NaN(size(fnum));
minp = NaN(size(fnum));

for i = 1:length(fnum)
    
    disp(i);
    
    p = pres(:,i);
    
    mxp = max(p); mnp = min(p);
    nn = ~isnan(p);
    
    p = p(nn);
    itop = find((mnp - lev) > 0,1,'last');
    if isempty(itop); itop = 1; end
    
    ibottom = find((mxp - lev) < 0,1,'first');
    if isempty(ibottom); ibottom = length(lev); end
    
    if (mnp - lev(itop)) > .1*(mnp + lev(itop))/2  || (mnp - lev(itop)) > 4
        itop = itop+1;
    end
    
    jjj = find(p < mxp,1,'last');
    
    if (mxp - lev(ibottom)) < -.2*(mxp + lev(ibottom))/2 || (mxp - lev(ibottom)) <= -28 || lev(ibottom)-mxp >= 2*(mxp - p(jjj))
        ibottom = ibottom-1;
    end
    
    maxp(i) = lev(ibottom);
    minp(i) = lev(itop);
    
end

clear i ibottom itop p i jjj mnp mxp nn p lev

D = dynht_cont(temp,pres,psal,level,reflevel,maxp,minp);

clear indnum maxp minp pres psal temp type reflevel

save /home/alison/Matlab/Mat' files'/Global' Mapping 2'/North' Pacific'/dynht_NP.mat

%% compute dynamic height values

clear; pack;

cd /home/alison/Matlab/M' files'/Global' Mapping 2'/

load /home/alison/Matlab/Mat' files'/Global' Mapping 2'/South' Pacific'/profiles_final.mat

reflevel = 900;
level = [5;10;20;30;50;75;100;125;150;200;250;300;400;500;600;700;800;900;1000;1100;1200;1300;1400;1500;1600;1700;1800;1900;2000];
lev = [0 5 10 20 30 50 75 100 125 150 175 200 250 275 300 350 400 450 500 550 600 650 700 750 800 850 900 950 1000 1050 1100 1150 1200 1250 1300 1350 1400 1450 1500 1550 1600 1650 1700 1750 1800 1850 1900 1950 2000]';

maxp = NaN(size(fnum));
minp = NaN(size(fnum));

for i = 1:length(fnum)
    
    disp(i);
    
    p = pres(:,i);
    
    mxp = max(p); mnp = min(p);
    nn = ~isnan(p);
    
    p = p(nn);
    itop = find((mnp - lev) > 0,1,'last');
    if isempty(itop); itop = 1; end
    
    ibottom = find((mxp - lev) < 0,1,'first');
    if isempty(ibottom); ibottom = length(lev); end
    
    if (mnp - lev(itop)) > .1*(mnp + lev(itop))/2  || (mnp - lev(itop)) > 4
        itop = itop+1;
    end
    
    jjj = find(p < mxp,1,'last');
    
    if (mxp - lev(ibottom)) < -.2*(mxp + lev(ibottom))/2 || (mxp - lev(ibottom)) <= -28 || lev(ibottom)-mxp >= 2*(mxp - p(jjj))
        ibottom = ibottom-1;
    end
    
    maxp(i) = lev(ibottom);
    minp(i) = lev(itop);
    
end

clear i ibottom itop p i jjj mnp mxp nn p lev

D = dynht_cont(temp,pres,psal,level,reflevel,maxp,minp);

clear indnum maxp minp pres psal temp type reflevel

save /home/alison/Matlab/Mat' files'/Global' Mapping 2'/South' Pacific'/dynht_SP.mat

%% compute dynamic height values

clear; pack;

cd /home/alison/Matlab/M' files'/Global' Mapping 2'/

load /home/alison/Matlab/Mat' files'/Global' Mapping 2'/Global/profiles_final_2004.mat

reflevel = 900;
level = [5;10;20;30;50;75;100;125;150;200;250;300;400;500;600;700;800;900;1000;1100;1200;1300;1400;1500;1600;1700;1800;1900;2000];
lev = [0 5 10 20 30 50 75 100 125 150 175 200 250 275 300 350 400 450 500 550 600 650 700 750 800 850 900 950 1000 1050 1100 1150 1200 1250 1300 1350 1400 1450 1500 1550 1600 1650 1700 1750 1800 1850 1900 1950 2000]';

maxp = NaN(size(fnum));
minp = NaN(size(fnum));

for i = 1:length(fnum)
    
    disp(i);
    
    p = pres(:,i);
    
    mxp = max(p); mnp = min(p);
    nn = ~isnan(p);
    
    p = p(nn);
    itop = find((mnp - lev) > 0,1,'last');
    if isempty(itop); itop = 1; end
    
    ibottom = find((mxp - lev) < 0,1,'first');
    if isempty(ibottom); ibottom = length(lev); end
    
    if (mnp - lev(itop)) > .1*(mnp + lev(itop))/2  || (mnp - lev(itop)) > 4
        itop = itop+1;
    end
    
    jjj = find(p < mxp,1,'last');
    
    if (mxp - lev(ibottom)) < -.2*(mxp + lev(ibottom))/2 || (mxp - lev(ibottom)) <= -28 || lev(ibottom)-mxp >= 2*(mxp - p(jjj))
        ibottom = ibottom-1;
    end
    
    maxp(i) = lev(ibottom);
    minp(i) = lev(itop);
    
end

clear i ibottom itop p i jjj mnp mxp nn p lev

D = dynht_cont(temp,pres,psal,level,reflevel,maxp,minp);

clear indnum maxp minp pres psal temp type reflevel

save /home/alison/Matlab/Mat' files'/Global' Mapping 2'/Global/dynht_2004.mat

