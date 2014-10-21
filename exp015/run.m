clear all
close all
% add library paths
restoredefaultpath
addpath(genpath('../../gsw_matlab_v3_02'))
addpath(genpath('../../omega/ansu_utils/external_scripts/'))
addpath(genpath('.'))

total_time=tic;
[s,ct,p,gamma_96,lats,longs]=gammanc_to_sctp;

[s,ct,gamma_96]=make_mask_96(s,ct,p,gamma_96,longs,lats);

save('data/input_data.mat')

user_input;


la=squeeze(lats(1,:,:));
lo=squeeze(longs(1,:,:));

[dy,dx]=scale_fac(la,lo);
%save('data/dxdy.mat', 'dx','dy') 
%load('data/dxdy.mat');
[nz,ny,nx]=size(s);
dx=repmat(permute(dx,[3 1 2]),[nz 1 1]);
dy=repmat(permute(dy,[3 1 2]),[nz 1 1]);
dz=circshift(p,[-1 0 0])-p;
dz(end,:,:)=nan;
dx=0*dx+1;
dy=0*dy+1;
dz=0*dz+1;
 
[ew99,ew99_mod]=ew99_modified(s,ct,p,dx,dy,dz);

display(['Total runtime ',num2str(toc(total_time)),' seconds'])


