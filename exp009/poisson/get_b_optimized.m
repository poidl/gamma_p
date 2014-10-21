function b=get_b_optimized(gp,n1,n2,n3,dx,dy,dz)

user_input;
[nz,ny,nx]=size(n1);

gam=~isnan(gp);

en= gam & circshift(gam,[ 0  0 -1]);
nn= gam & circshift(gam,[ 0 -1  0]);
ln= gam & circshift(gam,[-1  0  0]);

if ~zonally_periodic
    en(:,:,end)=false;
end
nn(:,end,:)=false;
ln(end,:,:)=false;

gam=gam(:);
% numbering well definied gammas
sreg=cumsum(gam);
sreg(~gam)=nan;

%ga=reshape(gam,[nz,ny,nx]);
sr=reshape(sreg,[nz,ny,nx]);

j_eg= circshift(sr,[0 0 -1]); % index of estern gridpoint
%j_wg= circshift(sr,[0 0  1]); % index of western gridpoint
j_ng= circshift(sr,[0 -1 0]); % north
%j_sg= circshift(sr,[0  1 0]); % south
%j_ug= circshift(sr,[ 1 0 0]); % upper
j_lg= circshift(sr,[-1 0 0]); % lower

j1_e=j_eg(en);
j2_e=sr(en);

j1_n=j_ng(nn);
j2_n=sr(nn);

j1_l=j_lg(ln);
j2_l=sr(ln);

i1e=1:length(j1_e); 
i1n=i1e(end)+(1:length(j1_n));
i1l=i1n(end)+(1:length(j1_l));

jcol=[j1_e;j1_n;j1_l;...
      j2_e;j2_n;j2_l];
   
irow=[i1e,i1n,i1l,...   
      i1e,i1n,i1l]';   
  
  
% n1e=circshift(n1,[0 0 -1]);
% n2n=circshift(n2,[0 -1 0]);
% n3l=circshift(n3,[-1 0 0]);

coef=0.5*[n1(en);n2(nn);n3(ln);...
          n1(en);n2(nn);n3(ln)];

A = sparse(irow,jcol,coef);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rhs
gpe=circshift(gp,[0 0 -1]);
gpn=circshift(gp,[0 -1 0]);
gpl=circshift(gp,[-1 0 0]);

gx=(gpe-gp);
gy=(gpn-gp);
gz=(gp-gpl);

y=[gx(en);gy(nn);gz(ln)];


n=sum(gam);
b_initial=ones(n,1);


nit=1000;

[b_,flag,relres,iter,resvec,lsvec] = lsqr(A,y,1e-15,nit,[],[],b_initial);

b=nan*gp;
b(gam)=b_;
%save_netcdf03(b,'b','data/b.nc')
%keyboard

