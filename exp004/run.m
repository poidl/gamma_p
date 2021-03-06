clear all
close all
% add library paths
restoredefaultpath
addpath(genpath('../../gsw_matlab_v3_02'))
addpath(genpath('../../omega/ansu_utils/external_scripts/'))
addpath(genpath('.'))

idealized01;
load('data/input_data_idealized.mat')



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

[divn,n1,n2,n3,interior,bdy,dx3,dy3,dz3]=div_n(s,ct,p,dx,dy,dz);
int=interior(:);
%sreg=cumsum(int); % label interior grid points
%sreg(~int)=nan;

% % boundary points which have a neighbouring intererior point
% ge= bdy & circshift(interior,[0  0  1]); % westward neighbour is interior pt -> form eastward gradient eq. here
% gw= bdy & circshift(interior,[0  0 -1]);
% gn= bdy & circshift(interior,[0  1  0]);
% gs= bdy & circshift(interior,[0 -1  0]);
% gl= bdy & circshift(interior,[1  0  0]);
% gu= bdy & circshift(interior,[-1 0  0]);
% boundary points which have a neighbour
ge= bdy & circshift(interior|bdy,[0  0  1]); 
gw= bdy & circshift(interior|bdy,[0  0 -1]);
gn= bdy & circshift(interior|bdy,[0  1  0]);
gs= bdy & circshift(interior|bdy,[0 -1  0]);
gl= bdy & circshift(interior|bdy,[1  0  0]);
gu= bdy & circshift(interior|bdy,[-1 0  0]);

if ~zonally_periodic
    ge(:,:,1)=false;
    gw(:,:,end)=false;
end
gn(:,1,:)=false;
gs(:,end,:)=false;
gl(1,:,:)=false;
gu(end,:,:)=false;

no_bdyeq= ge+gw+gn+gs+gl+gu; % number of boundary eq. at point
has_bdyeq=no_bdyeq(:)~=0;

gam= int | has_bdyeq;

% numbering well definied gammas
sreg=cumsum(gam);
sreg(~gam)=nan;

%ga=reshape(gam,[nz,ny,nx]);
sr=reshape(sreg,[nz,ny,nx]);

%upper=ga & circshift(ga,[1 0 0]);  upper(1,:)=false; % upper gridpoint exists
%lower=ga & circshift(ga,[-1 0 0]); lower(end,:)=false; % lower gridpoint exists


j_eg= circshift(sr,[0 0 -1]); % index of estern gridpoint
j_wg= circshift(sr,[0 0  1]); % index of western gridpoint
j_ng= circshift(sr,[0 -1 0]); % north
j_sg= circshift(sr,[0  1 0]); % south
j_ug= circshift(sr,[ 1 0 0]); % upper
j_lg= circshift(sr,[-1 0 0]); % lower


% TODO: THIS SHOULD BE A VOLUME INTEGRAL; WEIGHT BY VOLUME
dz3=0.5*(dz+circshift(dz,[1 0 0])); % re-grid onto original grid
dz3(1,:,:)=dz(1,:,:); % fix surface
for ii=1:nx*ny
    kk=find(isnan(dz3(:,ii)),1,'first');
    if ~isempty(kk) & kk~=1
        dz3(kk,ii)=dz(kk-1,ii);
    end
end
dy3(:,1,:)=dy(:,1,:); % fix south
dy3(:,end,:)=dy(:,end-1,:); % fix north
if ~zonally_periodic
    error('fix dz3')
end
vol=dx3.*dy3.*dz3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interior

j_e=j_eg(int);
j_w=j_wg(int);
j_n=j_ng(int);
j_s=j_sg(int);
j_u=j_ug(int);
j_l=j_lg(int);
j_c=sr(int);

dx_w=circshift(dx,[0 0 1]);
dy_s=circshift(dy,[0 1 0]);
dz_u=circshift(dz,[1 0 0]);

fe=vol(int)./(dx(int).^2); % x doesn't vary in zonal direction
fw=vol(int)./(dx_w(int).^2);
fn=vol(int)./(dy3(int).*dy_s(int));
fs=vol(int)./(dy3(int).*dy(int));
fu=vol(int)./(dz3(int).*dz_u(int));
fl=vol(int)./(dz3(int).*dz(int));

t1=(dy_s(int)+dy(int))./(dy_s(int).*dy(int).*dy3(int));
t2=(dz_u(int)+dz(int))./(dz_u(int).*dz(int).*dz3(int));
fc=-vol(int).*( 2./(dx3(int).^2) + t1+t2);

irow_int=bsxfun(@times,ones(1,7),(1:sum(int))');
irow_int=irow_int(:);
jcol_int=[j_e;j_w;j_n;j_s;j_u;j_l;j_c];

coef_int=[fe;fw;fn;fs;fu;fl;fc];
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% boundary
%keyboard

j1e=  sr(ge(:)); % plus one
j2e=j_wg(ge(:)); % minus one

j1w=j_eg(gw(:));
j2w=  sr(gw(:)); 

j1n=  sr(gn(:));
j2n=j_sg(gn(:));

j1s=j_ng(gs(:));
j2s=  sr(gs(:));

j1u=  sr(gu(:));
j2u=j_lg(gu(:));

j1l=j_ug(gl(:));
j2l=  sr(gl(:));


i1e=1:length(j1e); 
i1w=i1e(end)+(1:length(j1w));
i1n=i1w(end)+(1:length(j1n));
i1s=i1n(end)+(1:length(j1s));
i1u=i1s(end)+(1:length(j1u));
i1l=i1u(end)+(1:length(j1l));

jcol_bdy=[j1e;j1w;j1n;j1s;j1u;j1l; ...
          j2e;j2w;j2n;j2s;j2u;j2l];
 
irow_bdy=irow_int(end)+[i1e,i1w,i1n,i1s,i1u,i1l, ...   
                        i1e,i1w,i1n,i1s,i1u,i1l]';

dx_=dx(gam);
dy_=dy(gam);
dz_=dz(gam);
vol_=vol(gam);
coef_bdy=[ vol_(j2e)./dx_(j2e); vol_(j2w)./dx_(j2w); vol_(j2n)./dy_(j2n); vol_(j2s)./dy_(j2s); vol_(j1u)./dz_(j1u); vol_(j1l)./dz_(j1l); ... % ATTENTION: vert ax inverted
          -vol_(j2e)./dx_(j2e);-vol_(j2w)./dx_(j2w);-vol_(j2n)./dy_(j2n);-vol_(j2s)./dy_(j2s);-vol_(j1u)./dz_(j1u);-vol_(j1l)./dz_(j1l)];

n=sum(int)+sum(has_bdyeq);  

% condition
jcol_cond=(1:n)';     
irow_cond=(irow_bdy(end)+1)*ones(n,1);
coef_cond=(vol(gam)./sum(vol(gam))).*ones(n,1);
b_cond=0;

%jcol=jcol_int;
%irow=irow_int;
%coef=coef_int;
jcol=[jcol_int;jcol_bdy;jcol_cond];
irow=[irow_int;irow_bdy;irow_cond];
coef=[coef_int;coef_bdy;coef_cond];


A = sparse(irow,jcol,coef);

b_int=vol(int).*divn(int);
n1=vol(gam).*n1(gam);
n2=vol(gam).*n2(gam);
n3=vol(gam).*n3(gam);
b_bdy=[n1(j2e); n1(j2w); n2(j2n); n2(j2s); n3(j1u); n3(j1l)]; % ATTENTION: vert ax inverted 
b=[b_int;b_bdy;b_cond];
%b=b_int;

%keyboard
gamma_initial=zeros(n,1);
disp('starting LSQR()')
tic
[gamma,flag,relres,iter,resvec,lsvec] = lsqr(A,b,1e-15,1000,[],[],gamma_initial);
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

gamma_p=nan*s;
gamma_p(gam)=gamma;


keyboard
% no_eq= ge+gw+gn+gs+gl+gu; 
% has_eq=no_eq(:)~=0;
% no_eq=no_eq(has_eq); % number of equations at point
% nox=sum(has_eq);
% ibdy=sreg(has_eq);
% 
% jstart=1;
% c_e=1; % counter east
% c_w=1; % counter west
% c_n=1; % counter north
% c_s=1; % counter south
% 
% for ii=1:length(ibdy)
%     ix=ibdy(ii);
%     jend=jstart+no_eq(ix)-1;
%     j1(jstart:jend)=ix;
%     cnt=0;
%     if ge(ix)
%         i_e(c_e)=jstart+cnt;
%         cnt=cnt+1;
%         c_e=c_e+1;
%     end
%     if west(ix)
%         i_w(c_w)=jstart+cnt;
%         cnt=cnt+1;   
%         c_w=c_w+1;
%     end
%     if north(ix)
%         i_n(c_n)=jstart+cnt;
%         cnt=cnt+1;  
%         c_n=c_n+1;
%     end
%     if south(ix)
%         i_s(c_s)=jstart+cnt;
%         cnt=cnt+1; 
%         c_s=c_s+1;        
%     end
%     jstart=jend+1;    
% end
% %si=sum(int);
% 
% 
% keyboard
% 
% 
% 
