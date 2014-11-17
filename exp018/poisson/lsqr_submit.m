function gam=lsqr_submit(A,y,nit,gamma_initial)


disp('starting LSQR()')
tic
[gam,flag,relres,iter,resvec,lsvec] = lsqr(A,y,1e-15,nit,[],[],gamma_initial);
display(['LSQR() took ',num2str(toc),' seconds for ',num2str(length(lsvec)),' iterations']);
%keyboard
if length(lsvec)==length(resvec)
    mynorm=lsvec./resvec;
else
    mynorm=lsvec./resvec(2:end);
end
disp(['Arnorm/(anorm*rnorm) final: ', num2str(mynorm(end))])
disp(['Flag: ', num2str(flag)])
save('data/mynorm.mat','mynorm')

%save(['data/gamma_p_',num2str(ii),'.mat'])
