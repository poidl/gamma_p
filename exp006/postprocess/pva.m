figure()
h=imagesc(va)
set(h,'alphadata',~isnan(va))
set(gca,'ydir','normal')
colorbar()
