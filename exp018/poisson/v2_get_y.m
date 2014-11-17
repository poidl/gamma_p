function y=v2_get_y(n1,n2,n3,gam,j2e,j2w,j2n,j2s,j1u,j1l,b_cond)

%b_int=divn(int);
n1=n1(gam);
n2=n2(gam);
n3=n3(gam);
b_bdy=[n1(j2e); n1(j2w); n2(j2n); n2(j2s); n3(j1u); n3(j1l)]; % ATTENTION: vert ax inverted 
y=[b_bdy;b_cond];