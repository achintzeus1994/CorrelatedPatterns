clear all; close all; clc;
rng(3);
N     = 10000;
alpha = linspace(0.1,1,10);
K     = 50000*0.005;      % number of synapses to each neuron
p     = ceil(alpha*K);
c     = K/N;              % connection probability
len   = length(p);
m     = linspace(0,0.9,10);
len_m = length(m);
n_avg = 5;
ovlp_mem = zeros(len,len_m,n_avg);
for i = 1:len
    tic
    for j =1:len_m
        for k=1:n_avg
            ovlp_mem(i,j,k) = UlisesUncorrelated(i,p(i),N,m(j),k);
        end
    end
    toc
end
ovlp_mem1=zeros(len,len_m);
for i = 1:len
    for j =1:len_m
        ovlp_mem1(i,j) = mean(ovlp_mem(i,j,:));
    end
end
%% Calculating radius of attraction
m_val = zeros(len,1);
for i = 1 : len 
    cols  = find(ovlp_mem1(i,:)<0.6);
    if length(cols)>0
        m_val(i) = m(cols(1)); 
    elseif length(cols)==0
        m_val(i)=0.9;
    end
end
R = m_val;
%% Plotting figures
figure
subplot(2,1,1)
plot(alpha,ovlp_mem1(:,1),'-o','linewidth',2)
title('Overlap versus Load ($\rho=0$)','Interpreter','latex')
xlabel('Memory Load ($\alpha$)','Interpreter','latex');
ylabel('Overlap (m)','Interpreter','latex');
set(gca,'FontSize',18,'FontName','Times New Roman'); 
subplot(2,1,2)
plot(alpha,R,'-o','linewidth',2)
title('Radius of Attraction ($\rho=0$)','Interpreter','latex')
xlabel('Memory Load ($\alpha$)','Interpreter','latex');
ylabel('R','Interpreter','latex');
set(gca,'FontSize',18,'FontName','Times New Roman'); 
%%
save('Overlap_Matrix.mat','ovlp_mem')
%%
dat=load('Overlap_Matrix.mat');
%%
%corr_sim00 = [alpha',ovlp_mem];
%save('corr_sim00.mat','corr_sim00')