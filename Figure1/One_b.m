clear all; close all; clc;
Pc=2;
Pi=4;
rho = 0.5;      %correlation between patterns within each class
cormat_cat = rho^2*(ones(Pi-1)-eye(Pi-1))+eye(Pi-1);   %desired correlation matrix
f_row  = rho*ones(1,Pi-1);
f_col  = [1; rho*ones(Pi-1,1)];
cormat_cat = [f_row;cormat_cat];
cormat_cat = [f_col, cormat_cat];
cormat     = cormat_cat;
for i =1:Pc-1
    cormat = blkdiag(cormat,cormat_cat);
end
figure
imagesc(cormat)
title('Correlation Matrix ($\rho=0.5$)','Interpreter','latex'); % set title
colormap('jet'); % set the colorscheme
xlabel('Class 1                                   Class 2')
ylabel('Class 2                                   Class 1')
shading interp
set(gca,'FontSize',18,'FontName','Times New Roman');
axis square;
colormap(parula(10));
colorbar