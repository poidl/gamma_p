function b=get_b(gamma_p,n1,n2,n3)

user_input;
[nz,ny,nx]=size(n3);

n1=regrid(n1,3);
n2=regrid(n2,2);
n3=regrid(n3,1);


gpe=circshift(gamma_p,[0 0 -1]);
gpn=circshift(gamma_p,[0 -1 0]);
gpl=circshift(gamma_p,[-1 0 0]);

gx=(gpe-gamma_p);
gy=(gpn-gamma_p);
gz=(gamma_p-gpl);
% gx=(gpe-gamma_p)./dx;
% gy=(gpn-gamma_p)./dy;
% gz=(gamma_p-gpl)./dz;

if ~zonally_periodic
    gx(:,:,end)=nan;
end
gy(:,end,:)=nan;
gz(end,:,:)=nan;

gx=regrid(gx,3);
gy=regrid(gy,2);
gz=regrid(gz,1);

x=~isnan(n1) & ~isnan(gx);
y=~isnan(n2) & ~isnan(gy);
z=~isnan(n3) & ~isnan(gz);

b=nan*n1;
i1= x & ~y & ~z;
b(i1)=sqrt(gx(i1).^2)./sqrt(n1(i1).^2);

i1=~x &  y & ~z;
b(i1)=sqrt(gy(i1).^2)./sqrt(n2(i1).^2);

i1=~x & ~y &  z;
b(i1)=sqrt(gz(i1).^2)./sqrt(n3(i1).^2);

i1= x &  y & ~z;
b(i1)=sqrt(gx(i1).^2+gy(i1).^2)./sqrt(n1(i1).^2+n2(i1).^2);
i1= x & ~y &  z;
b(i1)=sqrt(gx(i1).^2+gz(i1).^2)./sqrt(n1(i1).^2+n3(i1).^2);
i1=~x &  y &  z;
b(i1)=sqrt(gy(i1).^2+gz(i1).^2)./sqrt(n2(i1).^2+n3(i1).^2);

i1=x & y & z;
b(i1)=sqrt(gx(i1).^2+gy(i1).^2+gz(i1).^2)./sqrt(n1(i1).^2+n2(i1).^2+n3(i1).^2);







% b3=gp_z./n3;
% 
% b1=b3;
% 
% for ii=1:nx*ny
%     kk=find(isnan(b1(:,ii)),1,'first');
%     if ~isempty(kk) & kk~=1
%         b1(kk,ii)=b1(kk-1,ii);
%     end
% end
% 
% b1=0.5*(b1+circshift(b1,[0 0 -1]));
% b2=0.5*(b1+circshift(b1,[0 -1 0]));
% if ~zonally_periodic
%     b1(:,:,end)=nan;
% end
% b2(:,end,:)=nan;
