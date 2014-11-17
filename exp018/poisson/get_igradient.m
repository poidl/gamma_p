function [ge,gw,gn,gs,gl,gu]=get_igradient(s,in)

user_input;

bdy=~isnan(s) & ~in;
ge= bdy & circshift(in, [ 0  0 -1]); % compute eastward gradient here (western bdy pt, eastern interior neighbour)
gw= bdy & circshift(in, [ 0  0  1]); 
gn= bdy & circshift(in, [ 0 -1  0]);
gs= bdy & circshift(in, [ 0  1  0]);
gl= bdy & circshift(in, [-1  0  0]);
gu= bdy & circshift(in, [ 1  0  0]);

if ~zonally_periodic
    ge(:,:,end)=false;
    gw(:,:,1)=false;
end
gn(:,end,:)=false;
gs(:,1,:)=false;
gl(end,:,:)=false;
gu(1,:,:)=false;