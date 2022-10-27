clear all; close all; clc;
rho    = 0.1;     %correlation between patterns within each class
Pc     = 1;     %number of classes
alow   = 0.1;   %low alpha
ahigh  = 0.8;   %high alpha
xi_ret = 1;     %retrived pattern
ntimes = ceil((ahigh-alow)*10+1);
N      = 10000;            % number of neurons           
K      = 50000*0.005;      % =250 from Ulises
c      = K/N;              % connection probability
alpha  = linspace(alow,ahigh,ntimes);
for i  = 1:ntimes
    Pi = floor(alpha*(c*N))+Pc;
end
ovlp   = zeros(ntimes,1);
for i  = 1:ntimes
    x       = Gutfreund_Rate(i,rho,Pc,Pi(i),xi_ret,N,c);
    x       = x(:,xi_ret);
    ovlp(i) = x(end);
end
%%
figure
plot(alpha,abs(ovlp),'-o','linewidth',2)
xlabel('Memory Load ($\alpha$)','Interpreter','latex');
ylabel('Overlap (m)','Interpreter','latex');
set(gca,'FontSize',18,'FontName','Times New Roman');
title('Simulation Analysis')
ylim([0 1]);
%% 
sIm_dat = [alpha;ovlp'];
save()
%%
%corr_sim00 = [alpha',ovlp_mem];
%save('corr_sim00.mat','corr_sim00')