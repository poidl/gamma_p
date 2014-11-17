function [ge,gw,gn,gs,gl,gu]=v2_get_igradient(s)

user_input;

nn=~isnan(s);
ge= nn & circshift(nn, [ 0  0 -1]); % compute eastward gradient here 
gw= nn & circshift(nn, [ 0  0  1]); 
gn= nn & circshift(nn, [ 0 -1  0]);
gs= nn & circshift(nn, [ 0  1  0]);
gl= nn & circshift(nn, [-1  0  0]);
gu= nn & circshift(nn, [ 1  0  0]);

if ~zonally_periodic
    ge(:,:,end)=false;
    gw(:,:,1)=false;
end
gn(:,end,:)=false;
gs(:,1,:)=false;
gl(end,:,:)=false;
gu(1,:,:)=false;
