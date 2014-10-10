a=[1:12]; a=reshape(a,[3,4]);

a=repmat(permute(a,[3 1 2]),[2 1 1]);
keyboard
ii= (a==4)|(a==8)|(a==6);

ii_=circshift(ii,[0 0 -1]);
a(ii)=a(ii_);

