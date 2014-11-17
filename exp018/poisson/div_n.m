function [divn,interior]=div_n(n1,n2,n3)

% south 
n2s=circshift(n2,[0 1 0]);
%dy3=0.5*(dy+circshift(dy,[0 1 0])); % re-grid onto original grid
% west
n1w=circshift(n1,[0 0 1]);
%dx3=0.5*(dx+circshift(dx,[0 0 1])); % re-grid onto original grid
% up
n3u=circshift(n3,[1 0 0]);
%dz3=0.5*(dz+circshift(dz,[1 0 0])); % re-grid onto original grid

dn1dx=(n1-n1w);
dn2dy=(n2-n2s);
dn3dz=(n3u-n3);

%dn1dx=(n1-n1w)./dx3;
%dn2dy=(n2-n2s)./dy3;
%dn3dz=(n3u-n3)./dz3;

interior= ~isnan(dn1dx) & ~isnan(dn2dy) & ~isnan(dn3dz);
divn=nan*n1;
divn(interior)=dn1dx(interior)+dn2dy(interior)+dn3dz(interior);



