clear all
close all
% add library paths
restoredefaultpath
addpath(genpath('../../gsw_matlab_v3_02'))
addpath(genpath('../../omega/ansu_utils/external_scripts/'))
addpath(genpath('.'))

idealized01;
load('data/input_data_idealized.mat')

total_time=tic;
%[s,ct,p,gamma_96,lats,longs]=gammanc_to_sctp;

%[s,ct,gamma_96]=make_mask_96(s,ct,p,gamma_96,longs,lats);

%save('data/input_data.mat')

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
 
[n1,n2,n3]=get_n(s,ct,p,dx,dy,dz);
[divn,x,y,z,xy,yz,xz,xyz,bdy_n,bdy_s,bdy_e,bdy_w,bdy_u,bdy_l,dx3,dy3,dz3]=div_n_new(n1,n2,n3,dx,dy,dz);
%bdy= ~isnan(s) & ~interior;
interior=x|y|z|xy|yz|xz|xyz;
int=interior(:);
%sreg=cumsum(int); % label interior grid points
%sreg(~int)=nan;

% boundary points 
ge= bdy_e; 
gw= bdy_w;
gn= bdy_n;
gs= bdy_s;
gl= bdy_l;
gu= bdy_u;


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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interior; 

dx_w=circshift(dx,[0 0 1]);
dy_s=circshift(dy,[0 1 0]);
dz_u=circshift(dz,[1 0 0]);

neqint=0;
%%%%%%% only in x direction
arg=x;

j_e=j_eg(arg);
j_w=j_wg(arg);
j_c=sr(arg);

fe=1./(dx(arg).^2); % x doesn't vary in zonal direction
fw=1./(dx_w(arg).^2);

fc=-(2./(dx3(arg).^2));

if any(arg(:))
    irow=bsxfun(@times,ones(1,3),(1:sum(arg(:)))');
    irow_int=irow(:);
    jcol_int=[j_e;j_w;j_c];

    coef_int=[fe;fw;fc];
else
    irow_int=[];
    jcol_int=[];
    coef_int=[];
end
neqint=neqint+sum(arg(:));

%%%%%%% only in y direction; 
arg=y;

j_n=j_ng(arg);
j_s=j_sg(arg);
j_c=sr(arg);

fn=1./(dy3(arg).*dy_s(arg));
fs=1./(dy3(arg).*dy(arg));

t1=(dy_s(arg)+dy(arg))./(dy_s(arg).*dy(arg).*dy3(arg));
fc=-t1;

neqint=0;
if any(arg(:))
    irow=bsxfun(@times,ones(1,3),(neqint+1:neqint+sum(arg(:)))');
    irow_int=[irow_int;irow(:)];
    jcol_int=[jcol_int;[j_n;j_s;j_c]];

    coef_int=[coef_int;[fn;fs;fc]];
else
    irow_int=[irow_int;[]];
    jcol_int=[jcol_int;[]];
    coef_int=[coef_int;[]];
end
neqint=neqint+sum(arg(:));

%%%%%%% only in z direction; 
arg=z;

j_u=j_ug(arg);
j_l=j_lg(arg);
j_c=sr(arg);

fu=1./(dz3(arg).*dz_u(arg));
fl=1./(dz3(arg).*dz(arg));

t2=(dz_u(arg)+dz(arg))./(dz_u(arg).*dz(arg).*dz3(arg));
fc=-t2;

if any(arg(:))
    irow=bsxfun(@times,ones(1,3),(neqint+1:neqint+sum(arg(:)))');
    irow_int=[irow_int;irow(:)];
    jcol_int=[jcol_int;[j_u;j_l;j_c]];

    coef_int=[coef_int;[fu;fl;fc]];
else
    irow_int=[irow_int;[]];
    jcol_int=[jcol_int;[]];
    coef_int=[coef_int;[]];
end
neqint=neqint+sum(arg(:));

%%%%%%% in x and y direction
arg=xy;

j_e=j_eg(arg);
j_w=j_wg(arg);
j_n=j_ng(arg);
j_s=j_sg(arg);
j_c=sr(arg);


fe=1./(dx(arg).^2); % x doesn't vary in zonal direction
fw=1./(dx_w(arg).^2);
fn=1./(dy3(arg).*dy_s(arg));
fs=1./(dy3(arg).*dy(arg));

t1=(dy_s(arg)+dy(arg))./(dy_s(arg).*dy(arg).*dy3(arg));
fc=-(2./(dx3(arg).^2)+t1);

if any(arg(:))
    irow=bsxfun(@times,ones(1,5),(neqint+1:neqint+sum(arg(:)))');
    irow_int=[irow_int;irow(:)];
    jcol_int=[jcol_int;[j_e;j_w;j_n;j_s;j_c]];

    coef_int=[coef_int;[fe;fw;fn;fs;fc]];
else
    irow_int=[irow_int;[]];
    jcol_int=[jcol_int;[]];
    coef_int=[coef_int;[]];
end
neqint=neqint+sum(arg(:));

%%%%%%% in y and z direction
arg=yz;

j_n=j_ng(arg);
j_s=j_sg(arg);
j_u=j_ug(arg);
j_l=j_lg(arg);
j_c=sr(arg);

fn=1./(dy3(arg).*dy_s(arg));
fs=1./(dy3(arg).*dy(arg));
fu=1./(dz3(arg).*dz_u(arg));
fl=1./(dz3(arg).*dz(arg));

t1=(dy_s(arg)+dy(arg))./(dy_s(arg).*dy(arg).*dy3(arg));
t2=(dz_u(arg)+dz(arg))./(dz_u(arg).*dz(arg).*dz3(arg));
fc=-(t1+t2);

if any(arg(:))
    irow=bsxfun(@times,ones(1,5),(neqint+1:neqint+sum(arg(:)))');
    irow_int=[irow_int;irow(:)];
    jcol_int=[jcol_int;[j_n;j_s;j_u;j_l;j_c]];

    coef_int=[coef_int;[fn;fs;fu;fl;fc]];
else
    irow_int=[irow_int;[]];
    jcol_int=[jcol_int;[]];
    coef_int=[coef_int;[]];
end
neqint=neqint+sum(arg(:));

%%%%%%% in x and z direction
arg=xz;

j_e=j_eg(arg);
j_w=j_wg(arg);
j_u=j_ug(arg);
j_l=j_lg(arg);
j_c=sr(arg);

fe=1./(dx(arg).^2); % x doesn't vary in zonal direction
fw=1./(dx_w(arg).^2);
fu=1./(dz3(arg).*dz_u(arg));
fl=1./(dz3(arg).*dz(arg));

t2=(dz_u(arg)+dz(arg))./(dz_u(arg).*dz(arg).*dz3(arg));
fc=-(2./(dx3(arg).^2)+t2);

if any(arg(:))
    irow=bsxfun(@times,ones(1,5),(neqint+1:neqint+sum(arg(:)))');
    irow_int=[irow_int;irow(:)];
    jcol_int=[jcol_int;[j_e;j_w;j_u;j_l;j_c]];

    coef_int=[coef_int;[fe;fw;fu;fl;fc]];
else
    irow_int=[irow_int;[]];
    jcol_int=[jcol_int;[]];
    coef_int=[coef_int;[]];
end
neqint=neqint+sum(arg(:));

%%%%%%% in all directions
arg=xyz;

j_e=j_eg(arg);
j_w=j_wg(arg);
j_n=j_ng(arg);
j_s=j_sg(arg);
j_u=j_ug(arg);
j_l=j_lg(arg);
j_c=sr(arg);

fe=1./(dx(arg).^2); % x doesn't vary in zonal direction
fw=1./(dx_w(arg).^2);
fn=1./(dy3(arg).*dy_s(arg));
fs=1./(dy3(arg).*dy(arg));
fu=1./(dz3(arg).*dz_u(arg));
fl=1./(dz3(arg).*dz(arg));

t1=(dy_s(arg)+dy(arg))./(dy_s(arg).*dy(arg).*dy3(arg));
t2=(dz_u(arg)+dz(arg))./(dz_u(arg).*dz(arg).*dz3(arg));
fc=-(2./(dx3(arg).^2)+t1+t2);

if any(arg(:))
    irow=bsxfun(@times,ones(1,7),(neqint+1:neqint+sum(arg(:)))');
    irow_int=[irow_int;irow(:)];
    jcol_int=[jcol_int;[j_e;j_w;j_n;j_s;j_u;j_l;j_c]];

    coef_int=[coef_int;[fe;fw;fn;fs;fu;fl;fc]];
else
    irow_int=[irow_int;[]];
    jcol_int=[jcol_int;[]];
    coef_int=[coef_int;[]];
end
neqint=neqint+sum(arg(:));

%keyboard


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interior;

% j_e=j_eg(ix);
% j_w=j_wg(ix);
% j_cx=sr(ix);
% 
% j_n=j_ng(iy);
% j_s=j_sg(iy);
% j_cy=sr(iy);
% 
% j_u=j_ug(iz);
% j_l=j_lg(iz);
% j_cx=sr(iz);
% 
% keyboard
% dx_w=circshift(dx,[0 0 1]);
% dy_s=circshift(dy,[0 1 0]);
% dz_u=circshift(dz,[1 0 0]);
% 
% fe=1./(dx(int).^2); % x doesn't vary in zonal direction
% fw=1./(dx_w(int).^2);
% fn=1./(dy3(int).*dy_s(int));
% fs=1./(dy3(int).*dy(int));
% fu=1./(dz3(int).*dz_u(int));
% fl=1./(dz3(int).*dz(int));
% 
% t1=(dy_s(int)+dy(int))./(dy_s(int).*dy(int).*dy3(int));
% t2=(dz_u(int)+dz(int))./(dz_u(int).*dz(int).*dz3(int));
% fc=-( 2./(dx3(int).^2) + t1+t2);
% 
% irow_int=bsxfun(@times,ones(1,7),(1:sum(int))');
% irow_int=irow_int(:);
% jcol_int=[j_e;j_w;j_n;j_s;j_u;j_l;j_c];
% 
% coef_int=[fe;fw;fn;fs;fu;fl;fc];
    

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
if ~isempty(i1e)
    i1w=i1e(end)+(1:length(j1w));
else
    i1w=1:length(j1w);
end
if ~isempty(i1w)
    i1n=i1w(end)+(1:length(j1n));
else
    i1n=1:length(j1n);
end
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

coef_bdy=[];
if ~isempty(j2e)
    coef_bdy=[ 1./dx_(j2e)];
end
if ~isempty(j2w)
    coef_bdy=[coef_bdy;1./dx_(j2w)];
end

coef_bdy=[coef_bdy; 1./dy_(j2n); 1./dy_(j2s); 1./dz_(j1u); 1./dz_(j1l)];
coef_bdy=[coef_bdy;-coef_bdy];
%coef_bdy=[ 1./dx_(j2e); 1./dx_(j2w); 1./dy_(j2n); 1./dy_(j2s); 1./dz_(j1u); 1./dz_(j1l); ... % ATTENTION: vert ax inverted
%          -1./dx_(j2e);-1./dx_(j2w);-1./dy_(j2n);-1./dy_(j2s);-1./dz_(j1u);-1./dz_(j1l)];

n=sum(gam);  

% condition
jcol_cond=(1:n)';     
irow_cond=(irow_bdy(end)+1)*ones(n,1);
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
dx3=dx3(gam);
dy3=dy3(gam);
dz3=dz3(gam);
vol=dx3.*dy3.*dz3;
coef_cond=(vol./sum(vol)).*ones(n,1);
b_cond=0;

%jcol=jcol_int;
%irow=irow_int;
%coef=coef_int;
jcol=[jcol_int;jcol_bdy;jcol_cond];
irow=[irow_int;irow_bdy;irow_cond];
coef=[coef_int;coef_bdy;coef_cond];


A = sparse(irow,jcol,coef);


yy=get_y(divn,n1,n2,n3,int,gam,j2e,j2w,j2n,j2s,j1u,j1l,b_cond);
%b=b_int;

%keyboard
gamma_initial=zeros(n,1);
nit=10000;
    %keyboard
    nit_p=20;
for ii=1:nit_p
    disp(['ii=',num2str(ii)]);
    disp('starting LSQR()')
    tic
    [gamma,flag,relres,iter,resvec,lsvec] = lsqr(A,yy,1e-15,nit,[],[],gamma_initial);
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
    save(['data/gamma_p_',num2str(ii),'.mat'])

    % lower
    bb=get_b(gamma_p,n1,n2,n3,dx,dy,dz);
        %keyboard
    n1=bb.*n1;
    n2=bb.*n2;
    n3=bb.*n3;
    [divn,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~]=div_n_new(n1,n2,n3,dx,dy,dz);
    yy=get_y(divn,n1,n2,n3,int,gam,j2e,j2w,j2n,j2s,j1u,j1l,b_cond);
    
    gamma_initial=gamma_p(gam);
    nit=400;
%     if ii==nit_p-1
%         nit=2000;
%     end
    if ii==1
        myb=bb;
    else
        myb=myb.*bb;
    end
    if ii==nit_p
        save_netcdf03(myb,'arg','data/b.nc')
        save_netcdf03(gamma_p,'gamma_p','data/gamma_p.nc')
    end
    
    %keyboard
end

display(['Total runtime ',num2str(toc(total_time)),' seconds'])


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
