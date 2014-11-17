function [irow_bdy,jcol_bdy,j2e,j2w,j2n,j2s,j1u,j1l]=matrix_ij_boundary(sr,ge,gw,gn,gs,gl,gu,...
                                       j_eg,j_wg,j_ng,j_sg,j_ug,j_lg)
j1e=j_eg(ge(:)); % plus one
j2e=  sr(ge(:)); % minus one

j1w=  sr(gw(:));
j2w=j_wg(gw(:)); 

j1n=j_ng(gn(:));
j2n=  sr(gn(:));

j1s=  sr(gs(:));
j2s=j_sg(gs(:));

j1u=j_ug(gu(:));
j2u=  sr(gu(:));

j1l=  sr(gl(:));
j2l=j_lg(gl(:));

i1e=1:length(j1e); 
if isempty(i1e) % possible in idealized, zonally periodic domains
    i1w=(1:length(j1w));
else
    i1w=i1e(end)+(1:length(j1w));
end
if isempty(i1e) & isempty(i1w) 
    i1n=(1:length(j1n));
elseif ~isempty(i1e) & ~isempty(i1w) 
    i1n=i1w(end)+(1:length(j1n));
else
    error('something''s wrong')
end
i1s=i1n(end)+(1:length(j1s));
i1u=i1s(end)+(1:length(j1u));
i1l=i1u(end)+(1:length(j1l));

jcol_bdy=[j1e;j1w;j1n;j1s;j1u;j1l; ...
          j2e;j2w;j2n;j2s;j2u;j2l];
 
irow_bdy=[i1e,i1w,i1n,i1s,i1u,i1l, ...   
         i1e,i1w,i1n,i1s,i1u,i1l]';