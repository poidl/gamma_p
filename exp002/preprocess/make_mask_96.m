function [s,ct,gam]=make_mask_96(s,ct,p,gam,lon,lat)

% gamma^n 96 mask
[nz,ny,nx]=size(s);
%keyboard
[ocean, n] = gamma_ocean_and_n(s,ct,p,lon,lat);

test1=ocean.*circshift(ocean,[-1 0]);
test2=ocean.*circshift(ocean,[0 -1]);

lo=lon(1,:)';
n= (test1(:)==5 | test2(:)==5) & (200<lo & lo<320);
setnan=n;

for kk=1:nz
    s(kk,setnan)=nan;
    ct(kk,setnan)=nan;
    gam(kk,setnan)=nan;
end


n= (test1==32);
setnan=n(:);

for kk=26:nz
    s(kk,setnan)=nan;
    ct(kk,setnan)=nan;
    gam(kk,setnan)=nan;
end



% h=imagesc(ocean)
% set(h,'alphadata',~isnan(ocean))
% set(gca,'ydir','normal')
% colorbar()
% 
% [ocean2, n2] = gamma_ocean_and_n(s,ct,p,lon,lat);
% figure()
% h=imagesc(ocean2)
% set(h,'alphadata',~isnan(ocean2))
% set(gca,'ydir','normal')
% colorbar()
% 
% figure
% va=squeeze(s(26,:,:));
% h=imagesc(va)
% set(h,'alphadata',~isnan(va))
% set(gca,'ydir','normal')
% colorbar()
