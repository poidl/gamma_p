function [b1,b2,b3]=get_b(gamma_p,n3,dz)
user_input;
[nz,ny,nx]=size(n3);

gpl=circshift(gamma_p,[-1 0 0]);
gp_z=(gamma_p-gpl)./dz;
gp_z(end,:,:)=nan;
b3=gp_z./n3;

b1=b3;

for ii=1:nx*ny
    kk=find(isnan(b1(:,ii)),1,'first');
    if ~isempty(kk) & kk~=1
        b1(kk,ii)=b1(kk-1,ii);
    end
end

b1=0.5*(b1+circshift(b1,[0 0 -1]));
b2=0.5*(b1+circshift(b1,[0 -1 0]));
if ~zonally_periodic
    b1(:,:,end)=nan;
end
b2(:,end,:)=nan;
