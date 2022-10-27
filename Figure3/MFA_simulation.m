%% Loading MFA data
clear all;close all;clc;
dat  = load('DataList.mat');
dat1 = load('DataListodd.mat');
dat2 = load('DataList50even.mat');
dat3 = load('DataList50odd.mat');
dat_rho00 = dat.Expression1;
dat_rho20 = dat.Expression2;
dat_rho40 = dat.Expression3;
dat_rho60 = dat.Expression4;

dat_rho10 = dat1.Expression1;
dat_rho30 = dat1.Expression2;
dat_rho50 = dat1.Expression3;
dat_rho70 = dat1.Expression4;

dat_rho52 = dat2.Expression1;
dat_rho54 = dat2.Expression2;
dat_rho56 = dat2.Expression3;
dat_rho58 = dat2.Expression4;

dat_rho51 = dat3.Expression1;
dat_rho53 = dat3.Expression2;
dat_rho55 = dat3.Expression3;
n_tot=16;
dat = zeros(9,n_tot);
dat(:,1) = dat_rho00(:,1);
dat(:,2) = dat_rho00(:,2);
dat(:,3) = dat_rho10(:,2);
dat(:,4) = dat_rho20(:,2);
dat(:,5) = dat_rho30(:,2);
dat(:,6) = dat_rho40(:,2);
dat(:,7) = dat_rho50(:,2);
dat(:,8) = dat_rho60(:,2);
dat(:,9) = dat_rho70(:,2);

dat(:,10) = dat_rho52(:,2);
dat(:,11) = dat_rho54(:,2);
dat(:,12) = dat_rho56(:,2);
dat(:,13) = dat_rho58(:,2);

dat(:,14) = dat_rho51(:,2);
dat(:,15) = dat_rho53(:,2);
dat(:,16) = dat_rho55(:,2);
figure
for i=2:9
    plot(dat(:,1),dat(:,i),'-o','linewidth',2)
    hold on
end
title('Mean Field Analysis')
axis tight
xlabel('Memory Load ($\alpha$)','Interpreter','latex');
ylabel('Overlap (m)','Interpreter','latex');
legend('$\rho$=0','$\rho$=0.1','$\rho$=0.2','$\rho$=0.3','$\rho$=0.4',...
    '$\rho$=0.5','$\rho$=0.6','$\rho$=0.7'...
    ,'Interpreter','latex')
set(gca,'FontSize',18,'FontName','Times New Roman'); 
% 
% figure
% plot(dat(:,1),dat(:,2),'-o','linewidth',2)
% hold on
% plot(sim00(:,1),sim00(:,2),'-o','linewidth',2)
% xlabel('Memory Load ($\alpha$)','Interpreter','latex');
% ylabel('Overlap (m)','Interpreter','latex');
% set(gca,'FontSize',18,'FontName','Times New Roman')