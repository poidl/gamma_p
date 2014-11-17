function [ew99]=ew99_modified(s,ct,p,lon,lat,modified)

user_input;

[dx,dy,dz]=get_dx_ones(lon,lat,p);
[nz,ny,nx]=size(s);

[n1,n2,n3]=get_n_old(s,ct,p,dx,dy,dz);
[divn,interior,dx3,dy3,dz3]=div_n_old(n1,n2,n3,dx,dy,dz);
bdy= ~isnan(s) & ~interior;
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

fe=1./(dx(int).^2); % x doesn't vary in zonal direction
fw=1./(dx_w(int).^2);
fn=1./(dy3(int).*dy_s(int));
fs=1./(dy3(int).*dy(int));
fu=1./(dz3(int).*dz_u(int));
fl=1./(dz3(int).*dz(int));

t1=(dy_s(int)+dy(int))./(dy_s(int).*dy(int).*dy3(int));
t2=(dz_u(int)+dz(int))./(dz_u(int).*dz(int).*dz3(int));
fc=-( 2./(dx3(int).^2) + t1+t2);

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
coef_bdy=[ 1./dx_(j2e); 1./dx_(j2w); 1./dy_(j2n); 1./dy_(j2s); 1./dz_(j1u); 1./dz_(j1l); ... % ATTENTION: vert ax inverted
          -1./dx_(j2e);-1./dx_(j2w);-1./dy_(j2n);-1./dy_(j2s);-1./dz_(j1u);-1./dz_(j1l)];

n=sum(int)+sum(has_bdyeq);  

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


b=get_y(divn,n1,n2,n3,int,gam,j2e,j2w,j2n,j2s,j1u,j1l,b_cond);
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
    [gamma,flag,relres,iter,resvec,lsvec] = lsqr(A,b,1e-15,nit,[],[],gamma_initial);
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
    bb=get_b_old(gamma_p,n1,n2,n3,dx,dy,dz);
    n1=bb.*n1;
    n2=bb.*n2;
    n3=bb.*n3;
    [divn,~,~,~,~]=div_n_old(n1,n2,n3,dx,dy,dz);
    b=get_y(divn,n1,n2,n3,int,gam,j2e,j2w,j2n,j2s,j1u,j1l,b_cond);
    
    gamma_initial=gamma_p(gam);
    nit=400;

    if ii==1
        myb=bb;
    else
        myb=myb.*bb;
    end
    if (ii==nit_p) || (ii==1 && modified==false)
        ew99=gamma_p;
        save_netcdf03(myb,'myb','data/bb.nc');
        save_netcdf03(gamma_p,'mygamma_p','data/gamma_p_.nc');
        break % for ii==1
    end

end
