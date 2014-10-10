close all
clear all

addpath(genpath('../../../../gsw_matlab_v3_02'))

xi=10;
yi=10;
zi=10;

ts=7; % surface temperature south
tn=5; % surface temperature north
dt=5; % surface temp minus bottom temp
dp=2e3;

y=linspace(0,1,yi);

tsn= ts+(tn-ts)*y;
tsn=repmat(tsn,[xi 1]);
tsn=repmat(permute(tsn,[3,2,1]),[zi 1 1]);

dt=-dt*linspace(0,1,zi);
dt=repmat(dt,[yi 1]);
dt=repmat(permute(dt,[2,1,3]),[1 1 xi]);

dp=dp*linspace(0,1,zi);
%dp=[0 200 1000];
dp=repmat(dp,[yi 1]);
p=repmat(permute(dp,[2,1,3]),[1 1 xi]);

ct=tsn+dt;
s=0*ct+35;


% h=contourf(squeeze(ct(:,:,1)))
% colorbar()
% set(gca,'YDir','reverse')
% figure()
% rpot=gsw_rho(s,ct,0*ct);
% contourf(squeeze(rpot(:,:,1)))
% colorbar()
% set(gca,'YDir','reverse')

%lats=[-84:4:84];
%lats=[-84:84:84];
lats=linspace(-84,84,yi);
lats=repmat(lats,[xi 1]);
lats=repmat(permute(lats,[3,2,1]),[zi 1 1]);

%longs=[0:4:356];
%longs=[0:178:356];
longs=linspace(0,356,xi);
longs=repmat(longs,[yi 1]);
longs=repmat(permute(longs,[3,1,2]),[zi 1 1]);

%keyboard
%set nans
%s(15)=nan; % 15 is center bottom
%ct(15)=nan; % 15 is center bottom
%s(24:25)=nan;
%ct(24:25)=nan;
%s(end-2:end,5,5)=nan;
%ct(end-2:end,5,5)=nan;
% s(3)=nan;
% ct(3)=nan;

vars = {'s','ct','p','lats','longs'};

save('data/input_data_idealized.mat',vars{:})





