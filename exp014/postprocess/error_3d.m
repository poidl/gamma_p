function [err_bar,err_med,values] = error_3d(va,sa,ct,p,vals)

    user_input;

    [nz,ny,nx]=size(va);

%     % find deepest location
%     pn=p; pn(isnan(sa))=nan;
%     imax=max(pn(:));
%     ibdy=int16(imax)/nz;
% 
%     bdy=va(:,ibdy);
%     bdy_top=bdy(1);
%     bdy_bottom=bdy(end);

    if nargin==4
        values=get_values(va);
    elseif nargin==5
        values=vals;
    end
           
    sbar=nan*ones(size(values));
    smed=nan*ones(size(values));

    for ii=1:length(values)
        [sx,sy,ss,cts,ps]=slope_error(va,sa,ct,p,values(ii));
        s=sx.^2+sy.^2;
        sbar(ii)=nanmean(s(:));
        smed(ii)=nanmedian(s(:));
    end
    
    err_bar=sbar;
    err_med=smed;
%     vap=sx.^2+sy.^2;
%     h=imagesc(vap);
%     set(h,'alphadata',~isnan(vap));
%     set(gca,'YDir','normal');
%     colorbar()
% 
% 
%     figure()
%     vap=ps
%     h=imagesc(vap);
%     set(h,'alphadata',~isnan(vap));
%     set(gca,'YDir','normal');
%     colorbar()
%     caxis([0 200])
end

function values=get_values(va)
    % get vector of values of va whose spacing decreases with
    % stratification (decreses with number of observations in a density range)

    va_min=min(va(:));
    va_max=max(va(:))+eps(1100); % add eps(1100) such that histc() returns zero for last bin.
    bins=linspace(va_min,va_max,30);
    hi=histc(va(:),bins);
    hi=hi(1:end-1); % remove last entry
    nsurf=100; % approx. number of total surfaces
    N=int16(hi*nsurf/sum(hi)); % number of surfaces for bin
    N=double(N);
    dg=diff(bins)./N'; % increment of gamma
    dg(~isfinite(dg))=nan;
    
    nsurf=sum(N);
    values=nan*(ones(1,nsurf));
    
    jj=1;
    for ii=1:length(dg);
        if ~isnan(dg(ii))
            values(jj: jj+N(ii)-1)=bins(ii)+dg(ii)*(0.5 : 1 : N(ii)-0.5);
            jj=jj+N(ii);
        end
    end 
end

function [sx,sy,ss,cts,ps]=slope_error(va,sa,ct,p,value)

    user_input;
    [nz,ny,nx]=size(sa);
    ps=var_on_surf_stef(p,va,value*ones(ny,nx));
    ss=var_on_surf_stef(sa,p,ps);
    cts=var_on_surf_stef(ct,p,ps);

    [ex,ey] = delta_tilde_rho(ss,cts,ps); % ex and ey are defined on staggered grids
    % regrid

    ex=0.5*(ex+circshift(ex,[0 1]));
    if ~zonally_periodic & nx~=1 % could extrapolate?
        ex(:,1)=nan;
        ex(:,end)=nan;
    end
    ey=0.5*(ey+circshift(ey,[1 0]));
    ey(1,:)=nan;
    ey(end,:)=nan;

    rho=gsw_rho(sa(:,:),ct(:,:),p(:,:));
    rho=reshape(rho,[nz,ny,nx]);
    rhos=var_on_surf_stef(rho,p,ps);
    [n2,pmid]=gsw_Nsquared(sa(:,:),ct(:,:),p(:,:));
    n2=reshape(n2,[nz-1,ny,nx]);
    pmid=reshape(pmid,[nz-1,ny,nx]);
    n2s=var_on_surf_stef(n2,pmid,ps);

    fac=(1/9.81)*rhos.*n2s;
    %facx=0.5*(fac+circshift(fac,[0 -1]));
    %facy=0.5*(fac+circshift(fac,[-1 0]));


    load('data/dxdy.mat')

    if nx~=1
        dx=0.5*(dx(:,1:end-1)+dx(:,2:end));
        dx=horzcat(dx(:,end), dx); % sloppy 
        ex=ex./dx;
    end
    dy=0.5*(dy(1:end-1,:)+dy(2:end,:));
    dy=vertcat(dy(end,:), dy); % sloppy
    
    ey=ey./dy;
    
    sx=ex./fac;
    sy=ey./fac;

end





















