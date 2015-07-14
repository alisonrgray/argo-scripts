%% south atlantic profile interpolation;

clear; pack

load /home/alison/Matlab/Ma't files'/Global' Mapping 2'/South' Atlantic'/prof_2005_qc.mat
fnum2 = fnum; latitude2 = latitude; longitude2 = longitude; time2 = time; type2 = type;
temp2 = NaN(996, 52815); psal2 = NaN(996, 52815); pres2 = NaN(996, 52815); datamode2 = datamode;
temp2_err = NaN(996, 52815); psal2_err = NaN(996, 52815); pres2_err = NaN(996, 52815);
[m,n] = size(pres);
temp2(1:m,1:n) = temp; psal2(1:m,1:n) = psal; pres2(1:m,1:n) = pres;
temp2_err(1:m,1:n) = temp_err; psal2_err(1:m,1:n) = psal_err; pres2_err(1:m,1:n) = pres_err;

load /home/alison/Matlab/M'at files'/Global' Mapping 2'/South' Atlantic'/prof_2006_qc.mat
fnum2 = [fnum2; fnum]; latitude2 = [latitude2; latitude]; longitude2 = [longitude2; longitude];
time2 = [time2; time]; type2 = [type2; type]; datamode2 = [datamode2; datamode];
[mm,nn] = size(pres);
temp2(1:mm,n+1:n+nn) = temp; psal2(1:mm,n+1:n+nn) = psal; pres2(1:mm,n+1:n+nn) = pres;
temp2_err(1:mm,n+1:n+nn) = temp_err; psal2_err(1:mm,n+1:n+nn) = psal_err; pres2_err(1:mm,n+1:n+nn) = pres_err;

load /home/alison/Matlab/M'at files'/Global' Mapping 2'/South' Atlantic'/prof_2007_qc.mat
fnum2 = [fnum2; fnum]; latitude2 = [latitude2; latitude]; longitude2 = [longitude2; longitude];
time2 = [time2; time]; type2 = [type2; type]; datamode2 = [datamode2; datamode];
[mmm,nnn] = size(pres);
temp2(1:mmm,n+nn+1:n+nn+nnn) = temp; psal2(1:mmm,n+nn+1:n+nn+nnn) = psal; pres2(1:mmm,n+nn+1:n+nn+nnn) = pres;
temp2_err(1:mmm,n+nn+1:n+nn+nnn) = temp_err; psal2_err(1:mmm,n+nn+1:n+nn+nnn) = psal_err; pres2_err(1:mmm,n+nn+1:n+nn+nnn) = pres_err;

load /home/alison/Matlab/M'at files'/Global' Mapping 2'/South' Atlantic'/prof_2008_qc.mat
fnum2 = [fnum2; fnum]; latitude2 = [latitude2; latitude]; longitude2 = [longitude2; longitude];
time2 = [time2; time]; type2 = [type2; type]; datamode2 = [datamode2; datamode];
[mmmm,nnnn] = size(pres);
temp2(1:mmmm,n+nn+nnn+1:n+nn+nnn+nnnn) = temp; psal2(1:mmmm,n+nn+nnn+1:n+nn+nnn+nnnn) = psal; pres2(1:mmmm,n+nn+nnn+1:n+nn+nnn+nnnn) = pres;
temp2_err(1:mmmm,n+nn+nnn+1:n+nn+nnn+nnnn) = temp_err; psal2_err(1:mmmm,n+nn+nnn+1:n+nn+nnn+nnnn) = psal_err; pres2_err(1:mmmm,n+nn+nnn+1:n+nn+nnn+nnnn) = pres_err;

load /home/alison/Matlab/M'at files'/Global' Mapping 2'/South' Atlantic'/prof_2009_qc.mat
fnum2 = [fnum2; fnum]; latitude2 = [latitude2; latitude]; longitude2 = [longitude2; longitude];
time2 = [time2; time]; type2 = [type2; type]; datamode2 = [datamode2; datamode];
[m5,n5] = size(pres);
temp2(1:m5,n+nn+nnn+nnnn+1:n+nn+nnn+nnnn+n5) = temp; psal2(1:m5,n+nn+nnn+nnnn+1:n+nn+nnn+nnnn+n5) = psal; pres2(1:m5,n+nn+nnn+nnnn+1:n+nn+nnn+nnnn+n5) = pres;
temp2_err(1:m5,n+nn+nnn+nnnn+1:n+nn+nnn+nnnn+n5) = temp_err; psal2_err(1:m5,n+nn+nnn+nnnn+1:n+nn+nnn+nnnn+n5) = psal_err; pres2_err(1:m5,n+nn+nnn+nnnn+1:n+nn+nnn+nnnn+n5) = pres_err;

load /home/alison/Matlab/M'at files'/Global' Mapping 2'/South' Atlantic'/prof_2010_qc.mat
fnum2 = [fnum2; fnum]; latitude2 = [latitude2; latitude]; longitude2 = [longitude2; longitude];
time2 = [time2; time]; type2 = [type2; type]; datamode2 = [datamode2; datamode];
[m6,n6] = size(pres);
temp2(1:m6,n+nn+nnn+nnnn+n5+1:n+nn+nnn+nnnn+n5+n6) = temp; psal2(1:m6,n+nn+nnn+nnnn+n5+1:n+nn+nnn+nnnn+n5+n6) = psal; pres2(1:m6,n+nn+nnn+nnnn+n5+1:n+nn+nnn+nnnn+n5+n6) = pres;
temp2_err(1:m6,n+nn+nnn+nnnn+n5+1:n+nn+nnn+nnnn+n5+n6) = temp_err; psal2_err(1:m6,n+nn+nnn+nnnn+n5+1:n+nn+nnn+nnnn+n5+n6) = psal_err; pres2_err(1:m6,n+nn+nnn+nnnn+n5+1:n+nn+nnn+nnnn+n5+n6) = pres_err;

clear fnum latitude longitude time type pres temp psal m mm mmm mmmm n nn nnn nnnn n5 n6 m5 m6

fnum = fnum2; latitude = latitude2; longitude = longitude2; type = type2; time = time2;
pres = pres2; psal = psal2; temp = temp2; datamode = datamode2;
pres_err = pres2_err; psal_err = psal2_err; temp_err = temp2_err;
clear fnum2 latitude2 longitude2 time2 type2 pres2 temp2 psal2 datamode2 pres2_err psal2_err temp2_err

clear errtype end_date param grayfloat start_date

ind = find(isnan(pres));
psal(ind) = NaN; temp(ind) = NaN;

ind = find(isnan(psal));
pres(ind) = NaN; temp(ind) = NaN;

ind = find(isnan(temp));
psal(ind) = NaN; pres(ind) = NaN;
clear ind

save /home/alison/Matlab/Mat' files'/Global' Mapping 2'/South' Atlantic'/profiles_interp.mat

%% remove double end values

k = 1; dblnum = [];
for i = 1:length(fnum)
    
    tmppres = pres(:,i);
    dp = tmppres(2:end)-tmppres(1:end-1);
    jj = find(~isnan(dp));
    
    if dp(jj(end)) < 10 && dp(jj(end-10)) > 10
        
        dblnum(k) = fnum(i); k = k+1; %#ok<SAGROW>
        
%         tmpsal = psal(:,i);
%         disp(fnum(i));
%         figure(1); clf; plot(tmpsal,-tmppres,'.-'); grid on; ylim([-pres(jj(end)+1,i)-50 -pres(jj(end)+1,i)+1000])
        
%         s = input('Enter 1 for bad, enter for good     ');
%         
%         if s     
%         pres(jj(end)+1,i) = NaN;
%         end
        %     else if dp(jj(end)) < 10 && dp(jj(end-10))<10
        %          tmpsal = psal(:,i);
        %         disp(fnum(i));
        %     figure(1); clf; plot(tmpsal,-tmppres,'.-'); grid on;
        % %     pause
        %         end
    end
    
    %     tmppres = pres(:,ind); tmpsal = psal(:,ind);
    %     figure(2); clf; plot(tmpsal(1:80,:),-tmppres(1:80,:),'.-'); grid on;
    
end

k = 1; indnum = [];
for i = 1:length(dblnum)
    
    if sum(dblnum(i) == indnum) == 0
        
     indnum(k) = dblnum(i); %#ok<SAGROW>
     k = k+1;
    end
end
indnum = indnum';
clear k i dp jj tmppres dblnum

for i = 1:length(indnum)
    
    
    ind = find(fnum == indnum(i));
    
%     figure(1); clf
%     plot(repmat((1:length(ind)),996,1),-pres(:,ind),'.-');
%     ylim([-max(max(pres(:,ind)))-200 -600]);
%     xlim([0 length(ind)+1]); grid on;
%     disp(indnum(i));
%     s = input('Enter 1 for remove double end values, enter for other    ');
    
%     if s
        for j = 1:length(ind)
            tmppres = pres(:,ind(j));
            dp = tmppres(2:end)-tmppres(1:end-1);
            jj = find(~isnan(dp));
            
            if dp(jj(end)) < 10 && dp(jj(end-10)) > 10
                pres(jj(end)+1,ind(j)) = NaN;
            end
            
            
        end
%     end
%     hold on;
%     plot(repmat((1:length(ind)),1015,1),-pres(:,ind),'m.-');
%     pause
%     if rem(i,10)==0
%         save /home/alison/Matlab/Mat' files'/Global' Mapping 2'/North' Atlantic'/tmp.mat
%     end
end

clear dp i ind j jj problemfloat tmppres tmpsal

ind = find(isnan(pres));
temp(ind) = NaN; psal(ind) = NaN;

clear ind indnum s

save /home/alison/Matlab/Mat' files'/Global' Mapping 2'/South' Atlantic'/profiles_interp.mat

%% examine for large gaps

lev = [0 5 10 20 30 50 75 100 125 150 175 200 250 275 300 350 400 450 500 550 600 650 700 750 800 850 900 950 1000 1050 1100 1150 1200 1250 1300 1350 1400 1450 1500 1550 1600 1650 1700 1750 1800 1850 1900 1950 2000]';

% 1) check for profiles with large gaps dp >= 400
indgap = []; k = 1;
for i = 1:length(fnum)
    p = pres(:,i);
    nn = ~isnan(p);
    p = p(nn);
    dp = p(2:end)-p(1:end-1);
    ind = find(dp >= 400, 1);
    if ~isempty(ind)
        indgap(k) = i; %#ok<SAGROW>
        k = k+1;
    end
end
indgap = indgap';
clear k i p nn dp ind

badind = NaN(length(indgap),1);
for i = 1:length(indgap)
    p = pres(:,indgap(i));
    nn = ~isnan(p);
    p = p(nn);
    dp = p(2:end)-p(1:end-1);
    ind = find(dp >= 500, 1);
    
    
    if p(ind) >= 1400  
        badind(i) = 0;
%         ii = find(pres(:,indgap(i)) == p(ind));
%         pres(ii+1:end,indgap(i)) = NaN;
    else if p(ind+1) <= 400
        badind(i) = 0;
        pause
        
        else
        badind(i) = 1;
        pres(:,indgap(i)) = NaN;
        end
    end
    
end

% view bad profiles;
% badind = NaN(length(indgap),1);
% ft = fittype('pchipinterp');
% 
% for j = 1:length(indgap)
%     
%     i = indgap(j);
%     p = pres(:,i);
%     
%     t = temp(:,i); s = psal(:,i);
%     mxp = max(p); mnp = min(p);
%     nn = ~isnan(p);
%     
%     p = p(nn);
%     dp = p(2:end)-p(1:end-1);
%     ind110 = find(dp >= 110);
%     jj = find( abs(time - time(i)) < 15 & fnum == fnum(i));
%     
%     t = t(nn);  s = s(nn);
%     
%     rho = sw_dens(s,t,p);
%     del = sw_svan(s,t,p);
%     
%     fxnt = fit(p,t,ft);
%     fxns = fit(p,s,ft);
%     fxnr = fit(p,rho,ft);
%     fxnd = fit(p,del,ft);
%     
%     for k = 1:length(ind110)
%                
%         itop = find((mnp - lev) > 0,1,'last');
%         if isempty(itop); itop = 1; end
%         %
%         ibottom = find((mxp - lev) < 0,1,'first');
%         if isempty(ibottom); ibottom = length(lev); end
%         
%         if abs(mnp - lev(itop)) > (mnp + lev(itop))/2*.1 || abs(mnp - lev(itop)) > 4
%             itop = itop+1;
%         end
%         
%         if abs(mxp - lev(ibottom)) > 30  %.5*(mxp + lev(ibottom))
%             ibottom = ibottom-1;
%         end
%         
%         rhojj = sw_dens(psal(:,jj),temp(:,jj),pres(:,jj));
%         deljj = sw_svan(psal(:,jj),temp(:,jj),pres(:,jj));
%         
%         figure(1); clf; plot(psal(:,jj),-pres(:,jj),'.-','color',[.8 .8 .8]);
%         hold on; plot(psal(:,i),-pres(:,i),'bo');
%         plot(feval(fxns,(lev(itop):5:lev(ibottom))),-(lev(itop):5:lev(ibottom)),'r-');
%         grid on;
%         
%         figure(2); clf; plot(temp(:,jj),-pres(:,jj),'.-','color',[.8 .8 .8]);
%         hold on; plot(temp(:,i),-pres(:,i),'bo');
%         plot(feval(fxnt,(lev(itop):5:lev(ibottom))),-(lev(itop):5:lev(ibottom)),'r-');
%         grid on;
%         
%         figure(3); clf; plot(rhojj,-pres(:,jj),'.-','color',[.8 .8 .8]);
%         hold on; plot(rho,-p,'bo');
%         plot(feval(fxnr,(lev(itop):5:lev(ibottom))),-(lev(itop):5:lev(ibottom)),'r-');
%         grid on;
%         
%         figure(4); clf; plot(deljj,-pres(:,jj),'.-','color',[.8 .8 .8]);
%         hold on; plot(del,-p,'bo');
%         plot(feval(fxnd,(lev(itop):5:lev(ibottom))),-(lev(itop):5:lev(ibottom)),'r-');
%         grid on;
%         
%         disp(fnum(i));
%         
%         %         if (sum(abs(mean(kkt,2)-kt)>1.5*std(kkt,[],2)) > 0 || min(size(kkt)) < 2)
%         %             disp('1');
%         %             dT(j,k) = 1;
%         %         else dT(j,k) = 0;
%         %             disp('0');
%         %         end
%         %         if (sum(abs(mean(kks,2)-ks)>1.5*std(kks,[],2)) > 0 || min(size(kks)) < 2)
%         %             disp('1');
%         %             dS(j,k) = 1;
%         %         else dS(j,k) = 0;
%         %             disp('0');
%         %         end
%     end
%     
%     problem = input('enter for good, 1 for bad');
%     if problem
%         badind(j) = 1;
%     else badind(j) = 0;
%     end
%     
% end

clear ans dp i i1 i2 ibottom ind110 ind2 itop j jj k kk kks kkt kp ks kt nn p problem tmps tmpt
clear dS dT badind del deljj ft fxnd fxnr fxns fxnt ii ind indgap mnp mxp rho rhoii s t rhojj

% remove profiles that are all NaNs
ind = find(sum(isnan(pres),1)~=996);

%%% removed 82 profiles for large gaps

fnum = fnum(ind); latitude = latitude(ind); longitude = longitude(ind);
pres = pres(:,ind); temp = temp(:,ind); psal = psal(:,ind); time = time(ind); type = type(ind);
datamode = datamode(ind); pres_err = pres_err(:,ind); temp_err = temp_err(:,ind); psal_err = psal_err(:,ind);
 
clear ind

save /home/alison/Matlab/Mat' files'/Global' Mapping 2'/South' Atlantic'/profiles_interp.mat

%% remove all profiles that don't include 900 db +/- 300 db

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

clear i ibottom itop p i jjj mnp mxp nn p 

% ind = find(~((minp <= 400 & maxp >= 900) | (minp <= 900 & maxp >= 1400)));
% for j = 1:length(ind)
%     
%     i = ind(j);
%     
%     figure(1); clf;
%     hold on; plot(psal(:,i),-pres(:,i),'bo');
% %     plot(fnval(fnxtr(fxn_pst{i,2}),(minp(i):5:maxp(i))),-(minp(i):5:maxp(i)),'r-');
%     grid on; axis tight;
%     
%     figure(2); clf; plot(temp(:,i),-pres(:,i),'bo'); hold on;
% %     plot(fnval(fnxtr(fxn_pst{i,1}),(minp(i):5:maxp(i))),-(minp(i):5:maxp(i)),'r-');
%     grid on; axis tight;
%     
%     disp(fnum(i))
% %     disp([minp(i) fxn_pst{i,3} maxp(i) fxn_pst{i,4}])
%     disp([latitude(i) longitude(i)])
%     
%     pause(.2)
%     %
%     %     problem = input('Enter for OK, 1 for fixable, 2 for not');
%     %
%     %     if problem == 1
%     %         badrho(j) = 1;
%     %     else if problem == 2
%     %             badrho(j) = 2;
%     %         else badrho(j) = 0;
%     %         end
%     %     end
%     %
%     
%     
% end
% 
% figure(1); clf
% scatter(longitude,latitude,'k');
% hold on;
% scatter(longitude(ind),latitude(ind),'r');

ind = find((minp <= 600 & maxp >= 900) | (minp <= 900 & maxp >= 1200));
latitude = latitude(ind); longitude = longitude(ind); time = time(ind);
fnum = fnum(ind); type = type(ind);
temp = temp(:,ind); pres = pres(:,ind); psal = psal(:,ind);
maxp = maxp(ind); minp = minp(ind); datamode = datamode(ind);
temp_err = temp_err(:,ind); pres_err = pres_err(:,ind); psal_err = psal_err(:,ind);

clear ind ans i j

save /home/alison/Matlab/Mat' files'/Global' Mapping 2'/South' Atlantic'/profiles_interp.mat


%% check for d(rho)/dp < 0

badind = NaN(size(fnum));
badpres = NaN(length(fnum),20);
ft = fittype('pchipinterp');
lev = [0 5 10 20 30 50 75 100 125 150 175 200 250 275 300 350 400 450 500 550 600 650 700 750 800 850 900 950 1000 1050 1100 1150 1200 1250 1300 1350 1400 1450 1500 1550 1600 1650 1700 1750 1800 1850 1900 1950 2000]';

for i = 1:length(fnum)
    p = pres(:,i);
    
    t = temp(:,i); s = psal(:,i);
    mxp = max(p); mnp = min(p);
    nn = ~isnan(p);
    
    p = p(nn);
    t = t(nn);  s = s(nn);
    
    rho = sw_dens(s,t,p);
    fxnr = fit(p,rho,ft);
    
    disp(i);
    
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
    
    tmp = differentiate(fxnr,lev(itop):1:lev(ibottom));
    
  
    if min(tmp) < -.02
        
        badind(i) = min(tmp);
        
        tmppres = find(tmp < -.02);
        
        badpres(i,1:length(tmppres)) = tmp(tmppres);
        
    end
end

clear fxnr i ibottom itop mnp mxp nn p s t rho tmppres tmp jjj 

indrho = find(~isnan(badind));
% badind = badind(indrho);
% badpres = badpres(indrho,:);

% look at profiles with drdp > 0.02
% for i = 1:length(indrho2)
%     j = indrho2(i);
%     
%     p = pres(:,j);
%     
%     t = temp(:,j); s = psal(:,j);
%     mxp = max(p); mnp = min(p);
%     nn = ~isnan(p);
%     
%     p = p(nn);
%     t = t(nn);  s = s(nn);
%     
%     rho = sw_dens(s,t,p);
%     fxnr = fit(p,rho,ft);
%     fxnt = fit(p,t,ft);
%     fxns = fit(p,s,ft);
%     
%     itop = find((mnp - lev) > 0,1,'last');
%     if isempty(itop); itop = 1; end
%     
%     ibottom = find((mxp - lev) < 0,1,'first');
%     if isempty(ibottom); ibottom = length(lev); end
%     
%     if (mnp - lev(itop)) > .1*(mnp + lev(itop))/2  || (mnp - lev(itop)) > 4
%         itop = itop+1;
%     end
%     
%     jjj = find(p < mxp,1,'last');
%     
%     if (mxp - lev(ibottom)) < -.2*(mxp + lev(ibottom))/2 || (mxp - lev(ibottom)) <= -28 || lev(ibottom)-mxp >= 2*(mxp - p(jjj))
%         ibottom = ibottom-1;
%     end
%     
%     figure(1); clf;
%     hold on; plot(psal(:,j),-pres(:,j),'bo');
%     plot(feval(fxns,(lev(itop):1:lev(ibottom))),-(lev(itop):1:lev(ibottom)),'r-');
%     grid on; axis tight;
%    
%     figure(2); clf; plot(temp(:,j),-pres(:,j),'bo'); hold on;
%     plot(feval(fxnt,(lev(itop):1:lev(ibottom))),-(lev(itop):1:lev(ibottom)),'r-');
%     grid on; axis tight;
%     
%     figure(3); clf;
%     plot(feval(fxnr,(lev(itop):1:lev(ibottom))),-(lev(itop):1:lev(ibottom)),'.-');
%     grid on; axis tight;
%     
%     
%     drdp = differentiate(fxnr,(lev(itop):1:lev(ibottom)));
%     figure(4); clf;
%     plot(drdp,-(lev(itop):1:lev(ibottom)),'b.-');
%     grid on; axis tight;
%     
%     disp(fnum(j));
%     problem = input('Enter for OK, 1 for fixable, 2 for not');
%     
%     if problem == 1
%         badrho(i) = 1; %#ok<SAGROW>
%     else if problem == 2
%             badrho(i) = 2; %#ok<SAGROW>
%         else badrho(i) = 0; %#ok<SAGROW>
%         end
%     end
%     
% end
% 
% indrho = indrho(badrho == 1);
% 
% clear badind badpres badrho drdp fxn* i ibottom itop j mnp mxp nn p s problem t rho tmp tmppres jjj

badind = NaN(length(indrho),2);

for i = 1:length(indrho)
    j = indrho(i);
    
    p = pres(:,j);
    
    t = temp(:,j); s = psal(:,j);
    mxp = max(p); mnp = min(p);
    nn = ~isnan(p);
    
    p = p(nn);
    t = t(nn);  s = s(nn);
    
    rho = sw_dens(s,t,p);
    fxnr = fit(p,rho,ft);
    fxnt = fit(p,t,ft);
    fxns = fit(p,s,ft);
    
    itop = find((mnp - lev) > 0,1,'last');
    if isempty(itop); itop = 1; end
    
    ibottom = find((mxp - lev) < 0,1,'first');
    if isempty(ibottom); ibottom = length(lev); end
    
    if (mnp - lev(itop)) > .1*(mnp + lev(itop))/2  || (mnp - lev(itop)) > 4
        itop = itop+1;
    end
    
    if (mxp - lev(ibottom)) < .2*(mxp + lev(ibottom))/2 || (mxp - lev(ibottom)) < -30
        ibottom = ibottom-1;
    end
    
    figure(1); clf;
    hold on; plot(psal(:,j),-pres(:,j),'bo');
    plot(feval(fxns,(lev(itop):1:lev(ibottom))),-(lev(itop):1:lev(ibottom)),'r-');
    grid on; axis tight;
   
    figure(2); clf; plot(temp(:,j),-pres(:,j),'bo'); hold on;
    plot(feval(fxnt,(lev(itop):1:lev(ibottom))),-(lev(itop):1:lev(ibottom)),'r-');
    grid on; axis tight;
    
    figure(3); clf;
    plot(feval(fxnr,(lev(itop):1:lev(ibottom))),-(lev(itop):1:lev(ibottom)),'.-');
    grid on; axis tight;
    
    drdp = differentiate(fxnr,(lev(itop):1:lev(ibottom)));
    figure(4); clf;
    plot(drdp,-(lev(itop):1:lev(ibottom)),'b.-');
    grid on; axis tight;
    
    disp(fnum(j));
    disp(pres(1:140,j)');
        
    problem = input('Enter values to delete   ');
    
    if ~isempty(problem)
    badind(i,1) = problem(1);
    badind(i,2) = problem(2);
    end
end

for i = 1:length(indrho)
   
    if ~isnan(badind(i,1))
       
        pres(badind(i,1):badind(i,2),indrho(i)) = NaN;
        
    end
    
end

ind = find(isnan(pres));
temp(ind) = NaN; psal(ind) = NaN;

clear ans badind badpres badrho drdp i ibottom ii ind indrho* itop j problem fxn*
clear nn mxp mnp p s t rho

% remove profiles that are all NaNs
ind = find(sum(isnan(pres),1)~=996);

fnum = fnum(ind); latitude = latitude(ind); longitude = longitude(ind);
pres = pres(:,ind); temp = temp(:,ind); psal = psal(:,ind); time = time(ind); type = type(ind);

clear ind

ind = find(isnan(pres));
temp(ind) = NaN; psal(ind) = NaN;

ind = find(isnan(temp));
pres(ind) = NaN; psal(ind) = NaN;

ind = find(isnan(psal));
pres(ind) = NaN; temp(ind) = NaN;

ind = find(isnan(pres));
temp_err(ind) = NaN; psal_err(ind) = NaN; pres_err(ind) = NaN;

clear ind ii jjj ft 

save /home/alison/Matlab/Mat' files'/Global' Mapping 2'/South' Atlantic'/profiles_interp.mat


%% find individual float numbers
%%%%%%%%%%%%%%%%%%
k = 1; indnum = 0;
for i = 1:length(fnum)
    if fnum(i) ~= indnum
        indnum(k) = fnum(i);
        k = k+1;
    end
end
indnum = indnum';
%%%%%%%%%%%%%%%%%%
clear i k


%% 
clear; pack;
load /home/alison/Matlab/Mat' files'/Global' Mapping 2'/South' Atlantic'/profiles_interp.mat
fnum2 = fnum; latitude2 = latitude; longitude2 = longitude; time2 = time; type2 = type;
temp2 = temp; psal2 = psal; pres2 = pres;

load /home/alison/Matlab/Mat' files'/Global' Mapping'/South' Atlantic'/profiles_interp.mat
fnum3 = fnum; latitude3 = latitude; longitude3 = longitude; time3 = time; type3 = type;
temp3 = temp; psal3 = psal; pres3 = pres;

load /home/alison/Matlab/Mat' files'/Global' Mapping'/South' Atlantic'/profiles.mat

clear float_num maxp minp indnum ft

ind = NaN(length(fnum2),1);
for i = 1:length(fnum2)
    p = pres2(:,i);
    ii = find(fnum3 == fnum2(i) & abs(latitude3 - latitude2(i)) < .5 & abs(longitude3 - longitude2(i)) < .5 & abs(time3-time2(i)) < .25);
    if length(ii) == 1
        pp = pres3(:,ii);
        if sum(~isnan(pp)) < sum(~isnan(p))
            ind(i) = 1;
        end
    else ind(i) = 0;
    end
end
clear i p ii pp 

ft = fittype('pchipinterp');
bad = NaN(length(fnum2),1);
for i = 50095:length(ind)
    if ind(i) == 1
        p = pres2(:,i); t = temp2(:,i); s = psal2(:,i);
        ii = find(fnum3 == fnum2(i) & abs(latitude3 - latitude2(i)) < .5 & abs(longitude3 - longitude2(i)) < .5 & abs(time3-time2(i)) < .25);
        pp = pres3(:,ii); tt = temp3(:,ii); ss = psal3(:,ii);
        
        d = sw_svan(s,t,p);
        dd = sw_svan(ss,tt,pp);
        
        figure(1); clf
        plot(d,-p,'.'); hold on;
        plot(dd,-pp,'r.');
        
        nn = find(~isnan(p));
        p = p(nn); d = d(nn);
        fxnd = fit(p,d,ft);
        plot(feval(fxnd,p(1):1:p(end)),-(p(1):1:p(end)),'k-');
        
        disp(max(p(2:end)-p(1:end-1)));
        
        b = input('enter for OK, 1 for delete same values as previously');
        if ~isempty(b)
            bad(i) = b;
        end
    end
end

for i = 1:length(bad)
    if bad(i) == 1
  
        p = pres2(:,i); t = temp2(:,i); s = psal2(:,i);
        ii = find(fnum3 == fnum2(i) & abs(latitude3 - latitude2(i)) < .5 & abs(longitude3 - longitude2(i)) < .5 & abs(time3-time2(i)) < .25);
        pp = pres3(:,ii); tt = temp3(:,ii); ss = psal3(:,ii);
        
%         d = sw_svan(s,t,p);
%         dd = sw_svan(ss,tt,pp);
%         
%         figure(1); clf
%         plot(d,-p,'.'); hold on;
%         plot(dd,-pp,'r.');
         
        iii = isnan(pp(1:996)) & ~isnan(p);
        pres2(iii,i) = NaN;
    end
    
end

clear ans b bad d dd ft fxnd i ii iii ind n nn p pp s ss t tt

ind = []; k = 1; l = 1; ind2 = []; ind3 = []; m = 1; ind4 = []; n = 1;
for i = 1:length(fnum2)
    ii = find(fnum3 == fnum2(i));
    iii = find(fnum == fnum2(i));
    
    if isempty(ii) && isempty(iii)
        ind3(m) = i; m = m+1;
    else if isempty(ii) && (type2(i) == 842 || (type2(i) == 847 || (type2(i) == 852 || type2(i) == 857)))
            ind4(n) = i; n = n+1;
        else if isempty(find(abs(latitude3(ii) - latitude2(i)) < .5 & abs(longitude3(ii) - longitude2(i)) < .5 & abs(time3(ii)-time2(i)) < .25, 1))
                
%                 disp(fnum2(i))
%                 figure(1); clf;
%                 plot(psal2(:,i),-pres2(:,i),'.-');
%                 figure(2); clf;
%                 plot(temp2(:,i),-pres2(:,i),'.-');
%                 pause
                
                iii = find(abs(latitude - latitude2(i)) < .5 & abs(longitude - longitude2(i)) < .5 & fnum == fnum2(i) & abs(time-time2(i)) < .25);
                if ~isempty(iii)
                    ind(k) = i;
                    k = k+1;
                else
                    ind2(l) = i;
                    l = l+1;
                    
                end
            end
        end
    end
end

clear l k i ii iii m n

%%% for ind3 = not in original dataset

k = 1; indnum = 0;
for i = 1:length(fnum2(ind3))
    if fnum2(ind3(i)) ~= indnum
        indnum(k) = fnum2(ind3(i));
        k = k+1;
    end
end
indnum = indnum';
%%%%%%%%%%%%%%%%%%
clear i k

for i = 1:length(indnum)
    ii = find(fnum2 == indnum(i));
    disp(indnum(i))
    figure(1); clf;
    plot(psal2(:,ii),-pres2(:,ii),'.-');
    figure(2); clf;
    plot(temp2(:,ii),-pres2(:,ii),'.-');
    pause
end

%%% for ind4 = fsi floats

k = 1; indnum = 0;
for i = 1:length(fnum2(ind4))
    if fnum2(ind4(i)) ~= indnum
        indnum(k) = fnum2(ind4(i));
        k = k+1;
    end
end
indnum = indnum';
%%%%%%%%%%%%%%%%%%
clear i k

for i = 1:length(indnum)
    ii = find(fnum2 == indnum(i));
    disp(indnum(i))
    figure(1); clf;
    plot(psal2(:,ii),-pres2(:,ii),'.-');
    figure(2); clf;
    plot(temp2(:,ii),-pres2(:,ii),'.-');
    pause
end

%%% for ind2 = profiles not in original dataset

k = 1; indnum = 0;
for i = 1:length(fnum2(ind2))
    if fnum2(ind2(i)) ~= indnum
        indnum(k) = fnum2(ind2(i));
        k = k+1;
    end
end
indnum = indnum';
%%%%%%%%%%%%%%%%%%
clear i k

for i = 1:length(indnum)
    ii = find(fnum2 == indnum(i));
    disp(indnum(i))
    figure(1); clf;
    plot(psal2(:,ii),-pres2(:,ii),'.-');
    figure(2); clf;
    plot(temp2(:,ii),-pres2(:,ii),'.-');
    pause
end

for i = 1:length(ind2)
    p = pres2(:,ind2(i));
    dp = p(2:end)-p(1:end-1);
    
    if max(dp) > 110
        figure(1); clf;
    plot(psal2(:,ind2(i)),-pres2(:,ind2(i)),'.-');
    figure(2); clf;
    plot(temp2(:,ind2(i)),-pres2(:,ind2(i)),'.-');
    disp(fnum2(ind2(i)))
    disp([max(dp) p(dp==max(dp))'])
    pause
    end
end


%%% for ind = deleted from original dataset

k = 1; indnum = 0;
for i = 1:length(fnum2(ind))
    if fnum2(ind(i)) ~= indnum
        indnum(k) = fnum2(ind(i));
        k = k+1;
    end
end
indnum = indnum';
%%%%%%%%%%%%%%%%%%
clear i k

for i = 1:length(indnum)
    ii = find(fnum2 == indnum(i));
    disp(indnum(i))
    figure(1); clf;
    plot(psal2(:,ii),-pres2(:,ii),'.-');
    figure(2); clf;
    plot(temp2(:,ii),-pres2(:,ii),'.-');
    pause
end

bad = NaN(length(ind),1);
for i = 1:length(ind)
    p = pres2(:,ind(i));
    dp = p(2:end)-p(1:end-1);
    
    if max(dp) > 110
        figure(1); clf;
    plot(psal2(:,ind(i)),-pres2(:,ind(i)),'.-');
    figure(2); clf;
    plot(temp2(:,ind(i)),-pres2(:,ind(i)),'.-');
    disp(fnum2(ind(i)))
    disp([max(dp) p(dp==max(dp))'])
    
    bad(i) = input('enter for good, 1 for bad');
    end
end

for i = 1:length(bad)
    if ~isnan(bad(i))
        figure(1); clf;
    plot(psal2(:,ind(i)),-pres2(:,ind(i)),'.-');
    figure(2); clf;
    plot(temp2(:,ind(i)),-pres2(:,ind(i)),'.-');
    disp(ind(i))
    pause
    end
end

pres2(:,347) = NaN;
pres2(:,696) = NaN;
pres2(:,1077) = NaN;
pres2(:,2817) = NaN;
pres2(1:3,2821) = NaN;
pres2(:,3556) = NaN;
pres2(:,9131) = NaN;
pres2(:,15072) = NaN;
pres2(:,33356) = NaN;
pres2(:,52073) = NaN;

clear pres pres3 temp temp3 fnum fnum3 i* dp bad ans latitude latitude3 longitude longitude3 p psal psal3 time time3 type type3
fnum = fnum2; latitude = latitude2; longitude = longitude2; pres = pres2; psal = psal2; temp = temp2; time = time2; type = type2;
clear pres2 temp2 psal2 time2 type2 latitude2 longitude2 fnum2

ind = find(isnan(pres));
temp(ind) = NaN; psal(ind) = NaN;

ind = find(isnan(temp));
pres(ind) = NaN; psal(ind) = NaN;

ind = find(isnan(psal));
pres(ind) = NaN; temp(ind) = NaN;

ind = find(isnan(pres));
temp_err(ind) = NaN; psal_err(ind) = NaN; pres_err(ind) = NaN;

ind = find(sum(isnan(pres),1)~=996);
fnum = fnum(ind); latitude = latitude(ind); longitude = longitude(ind);
pres = pres(:,ind); temp = temp(:,ind); psal = psal(:,ind); time = time(ind); type = type(ind);
datamode = datamode(ind); pres_err = pres_err(:,ind); temp_err = temp_err(:,ind); psal_err = psal_err(:,ind);

clear ind lev

save /home/alison/Matlab/Mat' files'/Global' Mapping 2'/South' Atlantic'/profiles_interp.mat

%%
ind = NaN(length(datamode),1);
for i = 1:length(datamode)
   
    if strcmp('A',datamode(i))
        ind(i) = 1;
    else if strcmp('R',datamode(i))
            ind(i) = 2;
        else
            ind(i) = 0;
        end
    end
end

Anum = fnum(ind==1);

k = 1; indnum = 0;
for i = 1:length(Anum)
    if Anum(i) ~= indnum
        indnum(k) = Anum(i);
        k = k+1;
    end
end
indnum = indnum';
%%%%%%%%%%%%%%%%%%
clear i k


for i = 1:length(indnum)
   
    ii = find(fnum == indnum(i));
    iii = find(fnum(logical(ind)) == indnum(i));
    disp(indnum(i))
    disp([length(ii) length(iii)]);
    indnum(i,2:3) = [length(ii) length(iii)];
    
    figure(1); clf
    plot(temp(:,ii),-pres(:,ii),'b.-');
    hold on;
%     plot(temp(:,iii),-pres(:,iii),'r.-');
    
    figure(2); clf
    plot(psal(:,ii),-pres(:,ii),'b.-');
    hold on;
%     plot(psal(:,iii),-pres(:,iii),'r.-');
%     pause
    
end

save /home/alison/Matlab/Mat' files'/Global' Mapping 2'/South' Atlantic'/profiles_final.mat

%% duplicates
tmp = [latitude longitude time];
dtmp = tmp(2:end,:) - tmp(1:end-1,:);
ind = find(dtmp(:,1) == 0 & dtmp(:,2) == 0 & dtmp(:,3) == 0);

i = 1;
pres(:,ind(i)) = NaN;

i = 2;
pres(:,ind(i)+1) = NaN;

i = 3;
pres(:,ind(i)) = NaN;

i = 4;
pres(:,ind(i)+1) = NaN;

i = 5;
pres(:,ind(i)+1) = NaN;

i = 6;
pres(:,ind(i)+1) = NaN;


fnum(ind(i))
t1 = temp(1:100,ind(i));
t2 = temp(1:100,ind(i)+1);
p1 = pres(1:100,ind(i));
p2 = pres(1:100,ind(i)+1);

latitude(ind(i))
 clf
plot(t1,-p1,'.-')
hold on;
plot(t2,-p2,'r.-');


ind = find(sum(isnan(pres),1)~=996);
fnum = fnum(ind); latitude = latitude(ind); longitude = longitude(ind);
pres = pres(:,ind); temp = temp(:,ind); psal = psal(:,ind); time = time(ind); type = type(ind);
datamode = datamode(ind); pres_err = pres_err(:,ind); temp_err = temp_err(:,ind); psal_err = psal_err(:,ind);

clear ind lev dtmp tmp t1 t2 p1 p2 i ans

save /home/alison/Matlab/Mat' files'/Global' Mapping 2'/South' Atlantic'/profiles_final.mat


