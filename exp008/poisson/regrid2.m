function a=regrid2(a,dim)


user_input;

[nz,ny,nx]=size(a);
%keyboard
if dim==3
elseif dim==2
    a=permute(a,[1,3,2]);
elseif dim==1
    a=permute(a,[2,3,1]);
end

aw=circshift(a,[0 0 1]);

anew=0.5*(aw+a); % onto original grid
if ~(dim==3 && zonally_periodic)
    anew(:,:,1)=a(:,:,1);
end

ii= isnan(anew) & ~isnan(a); % only consider eastern point here
anew(ii)=a(ii);

ii= isnan(anew) & ~isnan(aw); % only consider western point here
if ~(dim==3 && zonally_periodic)
    ii(:,:,end)=false;
end
anew(ii)=aw(ii);

a=anew;

if dim==2
    a=permute(a,[1,3,2]);
elseif dim==1
    a=permute(a,[3,1,2]);
end
