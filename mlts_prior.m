%To implement robust regression with least trimmed squares
% Parameters:
% x: P
% y: M
% P: P_rank
% M: M_rank
% alpha: timming proportion
function result=mlts_prior(x,y,P,M,alpha,options)

if nargin<6
    ns=1; nc=10; delta=0.01;
else
    ns=options.ns; nc=options.nc; delta=options.delta;
end

p=size(x,2); q=size(y,2); n=size(x,1);
h=floor(n*(1-alpha))+1;
obj0=1e10;

for i=1:ns
     xstart=P(2:end,:);ystart=M(2:end);
     bstart=pinv(xstart)*ystart;
    sigmastart=(ystart-xstart*bstart)'*(ystart-xstart*bstart)/q;
    for j=1:nc
        res=y-x*bstart;
        dist2=sum((res*pinv(sigmastart)).*res,2);
        [dist2s,idist2]=sort(dist2);
        idist2=idist2(1:h);
        xstart=x(idist2,:);ystart=y(idist2,:);
        bstart=pinv(xstart)*ystart;
        if h==p
          sigmastart=(ystart-xstart*bstart)'*(ystart-xstart*bstart)/((h-p)+0.01);
        else
           sigmastart=(ystart-xstart*bstart)'*(ystart-xstart*bstart)/(h-p);
        end
    end
    obj=det(sigmastart);
    if obj < obj0
        result.beta=bstart;
        result.sigma=sigmastart;
        obj0=obj;
    end
end

calpha = (1-alpha)/chi2cdf(chi2inv(1-alpha,q),q+2);
result.sigma=calpha*result.sigma;
res=y-x*result.beta;
result.dres=sum((res*pinv(result.sigma)).*res,2);
result.dres=(result.dres).^(0.5);


qdelta=sqrt(chi2inv(1-delta,q));
nooutlier=(result.dres<=qdelta*ones(n,1));
xgood=x(nooutlier,:);ygood=y(nooutlier,:);
result=pinv(xgood)*ygood;
