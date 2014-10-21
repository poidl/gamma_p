function [s,ct,p,gamma_96,lats,longs]=gammanc_to_sctp()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lat=ncread('/home/nfs/z3439823/mymatlab/gamma/gamma.nc','lat');
lon=ncread('/home/nfs/z3439823/mymatlab/gamma/gamma.nc','lon');
p=ncread('/home/nfs/z3439823/mymatlab/gamma/gamma.nc','pressure');
s=ncread('/home/nfs/z3439823/mymatlab/gamma/gamma.nc','s');
t=ncread('/home/nfs/z3439823/mymatlab/gamma/gamma.nc','t');
gamma=ncread('/home/nfs/z3439823/mymatlab/gamma/gamma.nc','gamma');

inan= s(:)==-99;
s(inan)=nan;
t(inan)=nan;
gamma(inan)=nan;

% lat=lat(2:end-1); % cut what is not in gk_ak_gamma.mat
% s=s(:,:,2:end-1);
% t=t(:,:,2:end-1);
% gamma=gamma(:,:,2:end-1);
% 
[nz nx ny]=size(s);
lat=repmat(lat,[1 nx]); lat=repmat(lat,[1 1 nz ]);
lon=repmat(lon,[1 ny]); lon=repmat(lon,[1 1 nz ]);
lats=permute(lat,[3 1 2]);
longs=permute(lon,[3 2 1]);

p=repmat(p,[1,nx]); p=repmat(p,[1 1 ny]);

s=permute(s, [1 3 2]);
t=permute(t, [1 3 2]);
p=permute(p, [1 3 2]);
gamma_96=permute(gamma, [1 3 2]);

%s=gsw_SA_from_SP(s,p,lon,lat);
ct=gsw_CT_from_t(s,t,p);

s=double(s);
ct=double(ct);
p=double(p);
gamma_96=double(gamma_96);

% wh=who;
% wh=setdiff(who,{'s','ct','p'});
% clear(wh{:});
% clear('wh');

end
