% dynamic height

% D = integral over pressure of specific volume.


function [D,P] = dynht_cont(temp,pres,psal,p0,reflevel,maxp,minp)

[~,m] = size(pres);
n = length(p0);

D = NaN(n,m);
P = NaN(n,m);

ft = fittype('pchipinterp');
 
for i = 1:m
    disp(i);
    p = pres(:,i); nn = find(~isnan(p)); p = p(nn);
    t = temp(:,i); t = t(nn);
    s = psal(:,i); s = s(nn);
    
    del = sw_svan(s,t,p);
    p = p*1e4;
    
    delspline = fit(p,del,ft);
    
    ind = find(p0 <= maxp(i) & p0 >= minp(i));
    pspline = p0(ind);
    
    dh = integrate(delspline,pspline*1e4,pspline(1)*1e4);
    
    dh = bsxfun(@plus,-dh,dh(pspline==reflevel));
    
    D(ind,i) = dh;
    P(ind,i) = pspline;
    
end

return

