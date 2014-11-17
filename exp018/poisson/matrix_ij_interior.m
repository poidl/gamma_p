function [irow_int,jcol_int]=matrix_ij_interior(int,sr,j_eg,j_wg,j_ng,j_sg,j_ug,j_lg)

j_e=j_eg(int);
j_w=j_wg(int);
j_n=j_ng(int);
j_s=j_sg(int);
j_u=j_ug(int);
j_l=j_lg(int);
j_c=sr(int);

irow_int=bsxfun(@times,ones(1,7),(1:sum(int))');
irow_int=irow_int(:);
jcol_int=[j_e;j_w;j_n;j_s;j_u;j_l;j_c];
