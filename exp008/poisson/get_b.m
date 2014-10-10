function [b,b1,b2,b3]=get_b(gamma_p,n1,n2,n3,dx,dy,dz)

    user_input;
    %[nz,ny,nx]=size(n3);

    gpe=circshift(gamma_p,[0 0 -1]);
    gpn=circshift(gamma_p,[0 -1 0]);
    gpl=circshift(gamma_p,[-1 0 0]);

    gx=(gpe-gamma_p)./dx;
    gy=(gpn-gamma_p)./dy;
    gz=(gamma_p-gpl)./dz;

    if ~zonally_periodic
        gx(:,:,end)=nan;
    end
    gy(:,end,:)=nan;
    gz(end,:,:)=nan;


    n1g=regrid2(n1,3);
    n2g=regrid2(n2,2);
    n3g=regrid2(n3,1);    
    
    gxg=regrid2(gx,3);
    gyg=regrid2(gy,2);
    gzg=regrid2(gz,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % interpolate onto n1 points

    n2t=itp(n2g,3);
    n3t=itp(n3g,3);
    
    gyt=itp(gyg,3);
    gzt=itp(gzg,3);   
    
    x=~isnan(gx);
    y=~isnan(gyt);
    z=~isnan(gzt);

    b1=b_func(gx,gyt,gzt,n1,n2t,n3t,x,y,z);

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % interpolate onto n2 points    
    
    n1t=itp(n1g,2);
    n3t=itp(n3g,2);
    
    gxt=itp(gxg,2);
    gzt=itp(gzg,2);   
    
    x=~isnan(gxt);
    y=~isnan(gy);
    z=~isnan(gzt);

    b2=b_func(gxt,gy,gzt,n1t,n2,n3t,x,y,z);

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % interpolate onto n3 points    
    
    n1t=itp(n1g,1);
    n2t=itp(n2g,1);
    
    gxt=itp(gxg,1);
    gyt=itp(gyg,1);   
    
    x=~isnan(gxt);
    y=~isnan(gyt);
    z=~isnan(gz);

    b3=b_func(gxt,gyt,gz,n1t,n2t,n3,x,y,z);   
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% forward differences

    n1=regrid(n1,3);
    n2=regrid(n2,2);
    n3=regrid(n3,1);

    gx=regrid(gx,3);
    gy=regrid(gy,2);
    gz=regrid(gz,1);

    x=~isnan(n1);
    y=~isnan(n2);
    z=~isnan(n3);

    b=b_func(gx,gy,gz,n1,n2,n3,x,y,z);

end

function b=b_func(gx,gy,gz,n1,n2,n3,x,y,z)

    b=nan*n1;
    i1= x & ~y & ~z;
    b(i1)=sqrt(gx(i1).^2)./sqrt(n1(i1).^2);

    i1=~x &  y & ~z;
    b(i1)=sqrt(gy(i1).^2)./sqrt(n2(i1).^2);

    i1=~x & ~y &  z;
    b(i1)=sqrt(gz(i1).^2)./sqrt(n3(i1).^2);

    i1= x &  y & ~z;
    b(i1)=sqrt(gx(i1).^2+gy(i1).^2)./sqrt(n1(i1).^2+n2(i1).^2);
    i1= x & ~y &  z;
    b(i1)=sqrt(gx(i1).^2+gz(i1).^2)./sqrt(n1(i1).^2+n3(i1).^2);
    i1=~x &  y &  z;
    b(i1)=sqrt(gy(i1).^2+gz(i1).^2)./sqrt(n2(i1).^2+n3(i1).^2);

    i1=x & y & z;
    b(i1)=sqrt(gx(i1).^2+gy(i1).^2+gz(i1).^2)./sqrt(n1(i1).^2+n2(i1).^2+n3(i1).^2);
    
    %b(b==Inf)=1; % necessary for Jackett/McDougall 97 data set.
end


function [v1]=itp(v1,dim)
    user_input;
    if dim==2
        v1=permute(v1,[1,3,2]);
    elseif dim==1
        v1=permute(v1,[2,3,1]);
    end

    v1=0.5*(v1+circshift(v1,[0 0 -1]));
    
    if ~(dim==3 && zonally_periodic)
        v1(:,:,end)=nan;
    end

    if dim==2
        v1=permute(v1,[1,3,2]);
    elseif dim==1
        v1=permute(v1,[3,1,2]);
    end
    
end
% b3=gp_z./n3;
% 
% b1=b3;
% 
% for ii=1:nx*ny
%     kk=find(isnan(b1(:,ii)),1,'first');
%     if ~isempty(kk) & kk~=1
%         b1(kk,ii)=b1(kk-1,ii);
%     end
% end
% 
% b1=0.5*(b1+circshift(b1,[0 0 -1]));
% b2=0.5*(b1+circshift(b1,[0 -1 0]));
% if ~zonally_periodic
%     b1(:,:,end)=nan;
% end
% b2(:,end,:)=nan;
