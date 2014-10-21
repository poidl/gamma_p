function a=regrid(a,dim)


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
in= isnan(a) & ~isnan(aw);
if ~(dim==3 && zonally_periodic)
    in(:,:,1)=false;
end
a(in)=aw(in);


if dim==2
    a=permute(a,[1,3,2]);
elseif dim==1
    a=permute(a,[3,1,2]);
end
