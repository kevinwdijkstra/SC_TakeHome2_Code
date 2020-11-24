%% Sparsity structure plot
n = 5;
[D2Mat_ssp, D3Mat_ssp] = CreateMatrix(n);

spy(D2Mat_ssp,'k.',8)
xlabel('')
set(gcf,'position',[100,100,400,400])
print -deps 2D_sparsity
close

spy(D3Mat_ssp,'k.',8)
xlabel('')
set(gcf,'position',[100,100,400,400])
print -deps 3D_sparsity
close

clear D2Mat_ssp D3Mat_ssp