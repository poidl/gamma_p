function [divn,n1,n2,n3,interior,boundary,dx3,dy3,dz3]=div_n(s,ct,p,dx,dy,dz)

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


% !! ASSUMING P VARIES IN Z ONLY !! %
n1=(gsw_rho(se,cte,p)-gsw_rho(s,ct,p));
%in=~isnan(s) & isnan(se);
%n1_=circshift(n1,[0 0 1]);
%n1(in)=n1_(in);
n1=n1./dx;

n2=(gsw_rho(sn,ctn,p)-gsw_rho(s,ct,p));
%n2(:,end,:)=n2(:,end-1,:);
%in=~isnan(s) & isnan(sn);
%n2_=circshift(n2,[0 1 0]);
%n2(in)=n2_(in);
n2=n2./dy;


pmid=0.5*(p+circshift(p,[-1 0 0]));
pmid(end,:,:)=nan;

n3=(gsw_rho(s,ct,pmid)-gsw_rho(sl,ctl,pmid));
%keyboard
n3=n3./dz;

if ~zonally_periodic
    n1(:,:,end)=nan;
end
n2(:,end,:)=nan;
n3(end,:,:)=nan;

if normalize
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
    n3=n3./(mag);
end

% south 
n2s=circshift(n2,[0 1 0]);
dy3=0.5*(dy+circshift(dy,[0 1 0])); % re-grid onto original grid
% west
n1w=circshift(n1,[0 0 1]);
dx3=0.5*(dx+circshift(dx,[0 0 1])); % re-grid onto original grid
% up
n3u=circshift(n3,[1 0 0]);
dz3=0.5*(dz+circshift(dz,[1 0 0])); % re-grid onto original grid

dn1dx=(n1-n1w)./dx3;
dn2dy=(n2-n2s)./dy3;
dn3dz=(n3u-n3)./dz3;

interior= ~isnan(dn1dx) & ~isnan(dn2dy) & ~isnan(dn3dz);
divn=nan*s;
divn(interior)=dn1dx(interior)+dn2dy(interior)+dn3dz(interior);

%boundary= (~isnan(dn1dx) | ~isnan(dn2dy) | ~isnan(dn3dz)) & ~interior;
boundary= ~isnan(s) & ~interior;


