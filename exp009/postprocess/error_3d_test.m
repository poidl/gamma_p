close all
clear all

addpath(genpath('../../../gsw_matlab_v3_02'))
addpath(genpath('..'))

% load SA CT p
load('../data/input_data.mat')
% load solution
load('../data/gamma_p_1.mat')

[s2bar,s2med,xdf]=error_3d(gamma_p,s,ct,p);
K=1e3;
dfbar=K*s2bar;
dfmed=K*s2med;

plot(xdf,dfbar,'r')
hold on
plot(xdf,dfbar,'ro')

plot(xdf,dfmed,'b')
plot(xdf,dfmed,'bo')
title('D_f [m^s/s]')
xlabel('\gamma^{rf} (black), \gamma^{p} (red)')

xlim([26.2,28])
%ylim([0 1.6e-8])
print('-dpng','-r200',['../figures/D_f.png'])

figure()
hist(gamma_p(:),50)
xlim([26.2,28])
save('../data/plots_error_3d_1.mat','xdf','dfbar','dfmed','gamma_p')

%print('-dpdf','-r200',['../figures/hist.pdf'])


