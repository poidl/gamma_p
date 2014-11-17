function [n1,n2,n3]=get_n(s,ct,p)

user_input;
normalize=false;

% north 
sn=circshift(s,[0 -1 0]);
ctn=circshift(ct,[0 -1 0]);

% east
se=circshift(s,[0 0 -1]);
cte=circshift(ct,[0 0 -1]);

% lower
sl=circshift(s,[-1 0 0]);
ctl=circshift(ct,[-1 0 0]);


% assuming p varies in z only
n1=(gsw_rho(se,cte,p)-gsw_rho(s,ct,p));
%in=~isnan(s) & isnan(se);
%n1_=circshift(n1,[0 0 1]);
%n1(in)=n1_(in);
%n1=n1./dx;

n2=(gsw_rho(sn,ctn,p)-gsw_rho(s,ct,p));
%n2(:,end,:)=n2(:,end-1,:);
%in=~isnan(s) & isnan(sn);
%n2_=circshift(n2,[0 1 0]);
%n2(in)=n2_(in);
%n2=n2./dy;


pmid=0.5*(p+circshift(p,[-1 0 0]));
pmid(end,:,:)=nan;

n3=(gsw_rho(s,ct,pmid)-gsw_rho(sl,ctl,pmid));
%n3=n3./dz;

if ~zonally_periodic
    n1(:,:,end)=nan;
end
n2(:,end,:)=nan;
n3(end,:,:)=nan;

if normalize % not really a normalization, rather use the zonal mean of n3 as vector magnitude
    [nz,ny,nx]=size(s);
    n3m=n3(1:end-1,:,:);
    n3m=repmat(nanmean(n3m,3),[1,1,nx]);  
    if min(n3m(:))>=-1e-8
        keyboard
        error('problem')
    end
    n3m=cat(1,n3m(1,:,:),n3m);
    mag=sqrt(n3m.^2);
    n1=n1./mag;
    n2=n2./mag;
    n3=n3./mag;

end





