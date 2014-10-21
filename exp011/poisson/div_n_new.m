function [divn,x,y,z,xy,yz,xz,xyz,bdy_n,bdy_s,bdy_e,bdy_w,bdy_u,bdy_l,dx3,dy3,dz3]=div_n_new(n1,n2,n3,dx,dy,dz)

user_input;

n2s=circshift(n2,[0 1 0]);
n1w=circshift(n1,[0 0 1]);
n3u=circshift(n3,[1 0 0]);

dx3=0.5*(dx+circshift(dx,[0 0 1])); % re-grid onto original grid
dy3=0.5*(dy+circshift(dy,[0 1 0])); 
dz3=0.5*(dz+circshift(dz,[1 0 0])); 

dn1dx=(n1-n1w)./dx3;
dn2dy=(n2-n2s)./dy3;
dn3dz=(n3u-n3)./dz3;

bdy_e= ( isnan(n1) & ~isnan(n1w));
bdy_w= (~isnan(n1) &  isnan(n1w));

bdy_n= ( isnan(n2) & ~isnan(n2s));
bdy_s= (~isnan(n2) &  isnan(n2s));

bdy_l= ( isnan(n3) & ~isnan(n3u));
bdy_u= (~isnan(n3) &  isnan(n3u));  
 
% THE FOLLOWING ARE NOT NECESSARY SINCE n1,n2,n3 are already set to nan at
% east/north/bottom

% if ~zonally_periodic
%     bdy_e(:,:,1)=false;
%     bdy_e(:,:,end)=true;
%     bdy_w(:,:,1)=true;  
%     bdy_w(:,:,end)=false;
% end
%bdy_n(:,1,:)=false;
%bdy_n(:,end,:)=true;
%bdy_s(:,1,:)=true;
%bdy_s(:,end,:)=false;
% bdy_l(1,:,:)=false;
% bdy_l(end,:,:)=true;
% bdy_u(1,:,:)=true;
% bdy_u(end,:,:)=false;
bdy=bdy_n|bdy_s|bdy_e|bdy_w|bdy_u|bdy_l;

%%%%%%%%%%%
ix= ~isnan(dn1dx) & ~bdy;
%keyboard
iy= ~isnan(dn2dy) & ~bdy;
iz= ~isnan(dn3dz) & ~bdy;
%keyboard
x=  ix & ~iy & ~iz;
y= ~ix &  iy & ~iz;
z= ~ix & ~iy &  iz;

xy=  ix &  iy & ~iz;
yz= ~ix &  iy &  iz;
xz=  ix & ~iy &  iz;

xyz= ix & iy & iz;

%%%%%%%%%%%

% interior= (~isnan(dn1dx) | ~isnan(dn2dy) | ~isnan(dn3dz)) & ~bdy;
% 
% x=interior & ( ~isnan(dn1dx) &  isnan(dn2dy) &  isnan(dn3dz) );
% y=interior & (  isnan(dn1dx) & ~isnan(dn2dy) &  isnan(dn3dz) );
% z=interior & (  isnan(dn1dx) &  isnan(dn2dy) & ~isnan(dn3dz) );
% 
% xy=interior & ( ~isnan(dn1dx) & ~isnan(dn2dy) &  isnan(dn3dz) );
% yz=interior & (  isnan(dn1dx) & ~isnan(dn2dy) & ~isnan(dn3dz) );
% xz=interior & ( ~isnan(dn1dx) &  isnan(dn2dy) & ~isnan(dn3dz) );
% 
% xyz=interior & ( ~isnan(dn1dx) & ~isnan(dn2dy) & ~isnan(dn3dz) );

divn=nan*n1;
divn(x)=dn1dx(x); 
divn(y)=dn2dy(y);
divn(z)=dn3dz(z);

divn(xy)=dn1dx(xy)+dn2dy(xy); 
divn(yz)=dn2dy(yz)+dn3dz(yz);
divn(xz)=dn1dx(xz)+dn3dz(xz);

divn(xyz)=dn1dx(xyz)+dn2dy(xyz)+dn3dz(xyz);




