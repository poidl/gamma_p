function [ew99]=v2_ew99_modified(s,ct,p,lon,lat,modified)
% don't solve Poisson equation, but solve min(\nabla\gamma-n) with LSQR.
% Note that this will yield a solution on all points that carry data, as
% opposed to the "clean" version which omits points which are not directly
% adjacent to the interior.

user_input;

[nz,ny,nx]=size(s);

[n1,n2,n3]=get_n(s,ct,p);
% in (interior) is true at points on which the full 3d Laplacian can be computed
% (6 neighbours)
%[divn,in]=div_n(n1,n2,n3);

[ge,gw,gn,gs,gl,gu]=v2_get_igradient(s); % compute eastward gradient here (has eastern neighbour)

%no_eq= ge+gw+gn+gs+gl+gu; % number of boundary eq. at point
gam=~isnan(s(:));

% numbering gammas
sreg=cumsum(gam);
sreg(~gam)=nan;

sr=reshape(sreg,[nz,ny,nx]);

j_eg= circshift(sr,[0 0 -1]); % index of estern gridpoint
j_wg= circshift(sr,[0 0  1]); % index of western gridpoint
j_ng= circshift(sr,[0 -1 0]); % north
j_sg= circshift(sr,[0  1 0]); % south
j_ug= circshift(sr,[ 1 0 0]); % upper
j_lg= circshift(sr,[-1 0 0]); % lower

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % interior
% [irow_int,jcol_int]=matrix_ij_interior(int,sr,j_eg,j_wg,j_ng,j_sg,j_ug,j_lg);
% 
% fe=ones(sum(int),1);
% fw=fe;fn=fe;fs=fe;fu=fe;fl=fe;
% fc=-6*fe;
% coef_int=[fe;fw;fn;fs;fu;fl;fc];
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% boundary
[irow_bdy,jcol_bdy,j2e,j2w,j2n,j2s,j1u,j1l]=matrix_ij_boundary(sr,ge,gw,gn,gs,gl,gu,...
                                       j_eg,j_wg,j_ng,j_sg,j_ug,j_lg);
%irow_bdy=irow_int(end)+irow_bdy;                          
%coef_bdy=[ 1./dx_(j2e); 1./dx_(j2w); 1./dy_(j2n); 1./dy_(j2s); 1./dz_(j1u); 1./dz_(j1l); ... % ATTENTION: vert ax inverted
%          -1./dx_(j2e);-1./dx_(j2w);-1./dy_(j2n);-1./dy_(j2s);-1./dz_(j1u);-1./dz_(j1l)];
coef_bdy=[ ge(ge); gw(gw); gn(gn); gs(gs); gu(gu); gl(gl); ... 
          -ge(ge);-gw(gw);-gn(gn);-gs(gs);-gu(gu);-gl(gl)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m=sum(gam);  % number of unknowns

% condition
jcol_cond=(1:m)';     
irow_cond=(irow_bdy(end)+1)*ones(m,1);
coef_cond=ones(m,1);
b_cond=0;

% jcol=jcol_int;
% irow=irow_int;
% coef=coef_int;
jcol=[jcol_bdy;jcol_cond];
irow=[irow_bdy;irow_cond];
coef=[coef_bdy;coef_cond];

A = sparse(irow,jcol,coef);

y=v2_get_y(n1,n2,n3,gam,j2e,j2w,j2n,j2s,j1u,j1l,b_cond);

gamma_initial=zeros(m,1);
nit=10000;
    %keyboard
    nit_p=20;
for ii=1:nit_p
    disp(['ii=',num2str(ii)]);
    gamma=lsqr_submit(A,y,nit,gamma_initial);    
    gamma_p=nan*s;
    gamma_p(gam)=gamma;

    b=get_b(gamma_p,n1,n2,n3);
    n1=b.*n1;
    n2=b.*n2;
    n3=b.*n3;
%     b=get_b_optimized(gamma_p,n1,n2,n3);
%     bx=0.5*(b+circshift(b,[0 0 -1]));
%     by=0.5*(b+circshift(b,[0 -1 0]));
%     bz=0.5*(b+circshift(b,[-1 0 0]));
%     n1=bx.*n1;
%     n2=by.*n2;
%     n3=bz.*n3;
%    [divn,~]=div_n(n1,n2,n3);
    y=v2_get_y(n1,n2,n3,gam,j2e,j2w,j2n,j2s,j1u,j1l,b_cond);
    
    gamma_initial=gamma_p(gam);
    nit=400;

    if ii==1
        myb=b;
    else
        myb=myb.*b;
    end
    if (ii==nit_p) || (ii==1 && modified==false)
        ew99=gamma_p;
        save_netcdf03(myb,'myb','data/b.nc');
        save_netcdf03(gamma_p,'gamma_p','data/gamma_p_.nc');
        break % for ii==1
    end

end
