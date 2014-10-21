close all
clear all

addpath(genpath('../../../gsw_matlab_v3_02'))
addpath(genpath('..'))

if 0
    % load SA CT p
    load('../data/input_data.mat')
    % load solution
    load('../data/gamma_p_20.mat')
else
    [s,ct,p,gamma_96,lats,longs]=gammanc_to_sctp;
    [s,ct,gamma_96]=make_mask_96(s,ct,p,gamma_96,longs,lats);
end

[s2bar,s2med,xdf]=error_3d(gamma_96,s,ct,p);
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
hist(gamma_96(:),50)
xlim([26.2,28])
save('../data/plots_error_3d_gamma_n.mat','xdf','dfbar','dfmed','gamma_96')

%print('-dpdf','-r200',['../figures/hist.pdf'])
