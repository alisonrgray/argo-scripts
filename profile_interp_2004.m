%% 2004 profile interpolation;

clear; pack

load /home/alison/Matlab/Mat' files'/Global' Mapping 2'/Global/prof_2004_qc.mat

ind = find(isnan(pres));
psal(ind) = NaN; temp(ind) = NaN;

ind = find(isnan(psal));
pres(ind) = NaN; temp(ind) = NaN;

ind = find(isnan(temp));
psal(ind) = NaN; pres(ind) = NaN;
clear ind

save /home/alison/Matlab/Mat' files'/Global' Mapping 2'/Global/profiles_interp_2004.mat

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
    ind = find(dp >= 400, 1);
    
    
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
ind = find(sum(isnan(pres),1)~=450);

%%% removed 82 profiles for large gaps

fnum = fnum(ind); latitude = latitude(ind); longitude = longitude(ind);
pres = pres(:,ind); temp = temp(:,ind); psal = psal(:,ind); time = time(ind); type = type(ind);
datamode = datamode(ind); pres_err = pres_err(:,ind); temp_err = temp_err(:,ind); psal_err = psal_err(:,ind);
 
clear ind

save /home/alison/Matlab/Mat' files'/Global' Mapping 2'/Global/profiles_interp_2004.mat

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

save /home/alison/Matlab/Mat' files'/Global' Mapping 2'/Global/profiles_interp_2004.mat

%% remove double end values

k = 1; dblnum = [];
for i = 1:length(fnum)
    
    tmppres = pres(:,i);
    dp = tmppres(2:end)-tmppres(1:end-1);
    jj = find(~isnan(dp));
    if length(jj) > 10
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
    else if dp(jj(end)) < 10 && dp(jj(end-8)) > 10
            
            dblnum(k) = fnum(i); k = k+1; %#ok<SAGROW,SAGRO
        end
        
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
    %     plot(repmat((1:length(ind)),1086,1),-pres(:,ind),'.-');
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
    %     plot(repmat((1:length(ind)),1086,1),-pres(:,ind),'m.-');
    %     pause
    %     if rem(i,10)==0
    %         save /home/alison/Matlab/Mat' files'/Global' Mapping'/North' Atlantic'/tmp.mat
    %     end
end

clear dp i ind j jj problemfloat tmppres tmpsal

ind = find(isnan(pres));
temp(ind) = NaN; psal(ind) = NaN;

clear ind indnum s

save /home/alison/Matlab/Mat' files'/Global' Mapping 2'/Global/profiles_interp_2004.mat

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
    
    if (mxp - lev(ibottom)) < -30  %.5*(mxp + lev(ibottom))
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
    problem = input('Enter for OK, 1 for fixable, 2 for not');
    
    if problem == 1
        badrho(i) = 1; %#ok<SAGROW>
    else if problem == 2
            badrho(i) = 2; %#ok<SAGROW>
        else badrho(i) = 0; %#ok<SAGROW>
        end
    end
    
end

indrho = indrho(badrho == 1);

clear badind badpres badrho drdp fxn* i ibottom itop j mnp mxp nn p s problem t rho tmp tmppres

clear ans badind badpres badrho drdp i ibottom ii ind indrho* itop j problem fxn*
clear nn mxp mnp p s t rho

% remove profiles that are all NaNs
ind = find(sum(isnan(pres),1)~=1001);

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

save /home/alison/Matlab/Mat' files'/Global' Mapping 2'/Global/profiles_interp_2004.mat



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

ind = find(ind == 1 | ind == 0);
fnum = fnum(ind); latitude = latitude(ind); longitude = longitude(ind);
pres = pres(:,ind); temp = temp(:,ind); psal = psal(:,ind); time = time(ind); type = type(ind);
datamode = datamode(ind); pres_err = pres_err(:,ind); temp_err = temp_err(:,ind); psal_err = psal_err(:,ind);



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

save /home/alison/Matlab/Mat' files'/Global' Mapping 2'/Global/profiles_final.mat

%% duplicates
tmp = [latitude longitude time fnum];
dtmp = tmp(2:end,:) - tmp(1:end-1,:);
ind = find(dtmp(:,1) == 0 & dtmp(:,2) == 0 & dtmp(:,3) == 0);

bad = NaN(20,1);
for i = 1:20
    t11 = t2(:,i);
    t21 = t3(:,i);
    s11 = s2(:,i);
    s21 = s3(:,i);
    p11 = p2(:,i);
    p21 = p3(:,i);
    figure(1); clf
    plot(t11,-p11,'.-') 
    hold on; plot(t21,-p21,'r.-');
    figure(2); clf
    plot(s11,-p11,'.-') 
    hold on; plot(s21,-p21,'r.-');
    
    disp(t11(1:5)');
    disp(t21(1:5)');
    disp([d2(i) d3(i)]);
    bad(i) = input('1 for delete blue, 2 for delete red');
end

for i = 1:length(ind)
    pres(:,ind(i)) = NaN;
end

ind = find(sum(isnan(pres),1)~=1001);
fnum = fnum(ind); latitude = latitude(ind); longitude = longitude(ind);
pres = pres(:,ind); temp = temp(:,ind); psal = psal(:,ind); time = time(ind); type = type(ind);
datamode = datamode(ind); pres_err = pres_err(:,ind); temp_err = temp_err(:,ind); psal_err = psal_err(:,ind);

clear ind lev dtmp tmp t1 t2 p1 p2 i ans bad s2 s1

[~,ind] = unique([latitude longitude time],'rows');
fnum = fnum(ind); latitude = latitude(ind); longitude = longitude(ind);
pres = pres(:,ind); temp = temp(:,ind); psal = psal(:,ind); time = time(ind); type = type(ind);
datamode = datamode(ind); pres_err = pres_err(:,ind); temp_err = temp_err(:,ind); psal_err = psal_err(:,ind);

clear ind

save /home/alison/Matlab/Mat' files'/Global' Mapping 2'/Global/profiles_final_2004.mat

%%
for i = 1:length(fnum)
p = pres(:,i);
p = p(~isnan(p));
maxdp(i) = max(p(2:end)-p(1:end-1));
end
ind = find(maxdp >= 110);

bad = NaN(length(ind),1);
for i = 22:length(ind)
figure(1); clf
plot(temp(:,ind(i)),-pres(:,ind(i)),'.-');
grid on;
figure(2); clf
plot(psal(:,ind(i)),-pres(:,ind(i)),'.-');
grid on;
b = input('bad? enter for OK, 1 for bad');
if ~isempty(b)
bad(i) = 1;
end
end

ind = ind(bad == 1);
clear bad
for i = 7:length(ind)
figure(1); clf
plot(temp(:,ind(i)),-pres(:,ind(i)),'.-');
grid on;
figure(2); clf
plot(psal(:,ind(i)),-pres(:,ind(i)),'.-');
grid on;
b = input('bad pressure range? ');
bad(i,:) = b;
end

for i = 1:length(ind)
   p = pres(:,ind(i));
   ii = find(p <= bad(i,2) & p >= bad(i,1));
   pres(ii,ind(i)) = NaN;
end

clear ans b bad i ii ind maxdp p

% remove profiles that are all NaNs
ind = find(sum(isnan(pres),1)~=450);

fnum = fnum(ind); latitude = latitude(ind); longitude = longitude(ind);
pres = pres(:,ind); temp = temp(:,ind); psal = psal(:,ind); time = time(ind); type = type(ind);
datamode = datamode(ind); pres_err = pres_err(:,ind); temp_err = temp_err(:,ind); psal_err = psal_err(:,ind);

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

ind = find(~isnan(latitude));
fnum = fnum(ind); latitude = latitude(ind); longitude = longitude(ind);
pres = pres(:,ind); temp = temp(:,ind); psal = psal(:,ind); time = time(ind); type = type(ind);
datamode = datamode(ind); pres_err = pres_err(:,ind); temp_err = temp_err(:,ind); psal_err = psal_err(:,ind);

clear ind

save /home/alison/Matlab/Mat' files'/Global' Mapping 2'/Global/profiles_final_2004.mat


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




