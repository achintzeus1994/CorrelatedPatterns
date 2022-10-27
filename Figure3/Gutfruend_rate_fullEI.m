clear all; close all; clc;
%% Parameter initialization
t_E  = 20;               % in ms 
t_I  = 5;
N_E  = 10000;           % number of excitatory neurons
N_I  = ceil(N_E/4);     % number of inhibitory neurons
c_EE = 2*0.025;         % EE connection probability
c_EI = 0.02;
c_IE = 0.02;
K_EE = N_E*c_EE;
K_EI = N_I*c_EI;
K_IE = N_E*c_IE;
tspan= 1000;             % simulation time= 2.5 sec
dt   = 0.5;              % time step 0.5 ms
A_EE = 1.5*4;
time = dt*linspace(0,tspan/dt,tspan/dt+1);
Tlen = length(time);
t_stim_on  = 0.4*Tlen;
t_stim_off = 0.6*Tlen;
Pc  = 5;        %number of classes
Pi  = 5;        %number of individual patterns/class
p   = Pc*Pi;
rho = 0.25;      %correlation between patterns within each class
cormat_cat = rho^2*(ones(Pi-1)-eye(Pi-1))+eye(Pi-1);   %desired correlation matrix
f_row  = rho*ones(1,Pi-1);
f_col  = [1; rho*ones(Pi-1,1)];
cormat_cat = [f_row;cormat_cat];
cormat_cat = [f_col, cormat_cat];
cormat     = cormat_cat;
for i =1:Pc-1
    cormat = blkdiag(cormat,cormat_cat);
end
%% Sample neurons
nsam = 20;
rneuE = randi(N_E,nsam,1);
rsamE = zeros(nsam,Tlen);
rneuI = randi(N_I,nsam,1);
rsamI = zeros(nsam,Tlen);
%% Recurrent matrix
nsypEE = zeros(N_E,K_EE);
nsypEI = zeros(N_E,K_EI);
nsypIE = zeros(N_I,K_IE);
for i =1:N_E
    nsypEE(i,:) = randperm(N_E,K_EE);    % number of input synapses/neurons
    nsypEI(i,:) = randperm(N_I,K_EI);
end
for i = 1:N_I
    nsypIE(i,:) = randperm(N_E,K_IE);
end
nneuEE  = repmat(linspace(1,N_E,N_E)',1,K_EE);    
nneuEI  = repmat(linspace(1,N_E,N_E)',1,K_EI);
nneuIE  = repmat(linspace(1,N_I,N_I)',1,K_IE);
cm_EE   = sparse(nneuEE,nsypEE,1,N_E,N_E);
cm_EI   = sparse(nneuEI,nsypEI,1,N_E,N_I);
cm_IE   = sparse(nneuIE,nsypIE,1,N_I,N_E);
clear vars nneuEE nneuEI nneuIE nsypEE nsypEI nsypIE
%% Loading pattern
L  = chol(cormat);                         %Cholesky decomposition of correlation matrix
xi = null([rand(N_E-(p+1),N_E);ones(1,N_E)])*L;
xi = xi - mean(xi);
xi = xi./std(xi);
save('xi_rho_25_5pc_5pi.mat','xi'); 
%xi = load('xi_rho_30_5pc_5pi.mat'); xi = xi.xi;
y= corrcoef(xi); y(y<0.00001)=0;
xic=xi(:,1);
for i=Pi:Pi:Pc*Pi-1
    xic=[xic,xi(:,i+1)];
end
del_var = 0:Pi:Pc*Pi-1;
xi(:,del_var+1)=[];
%% Plotting pattern correlations
rho_p = zeros(p-Pc,1);
flg   = 1;
for i=1:Pc
    for j =1:Pi-1
        rho_p(flg) = corr(f(sig(xic(:,i))),f(sig(xi(:,j+(i-1)*(Pi-1)))));
        flg = flg + 1;
    end
end
rho_p = mean(rho_p);
clear vars flg
%% Creating EE weight matrix
thrd     = 0.12;
xidiff_f = f(sig(xi))-rho_p*f(sig(repelem(xic,1,Pi-1)));
xidiff_g = g(sig(xi))-rho_p*g(sig(repelem(xic,1,Pi-1)));
J_EE = cm_EE/sqrt(K_EE).*max(A_EE/sqrt(K_EE)*xidiff_f*xidiff_g'+thrd,0);
for i = 1:N_E
    J_EE(i,i)=0;  % Getting rid of autapses
end
%% Calculating EI condition
J_IE  = 0.0673;
J_EI  = 0.0250;
alpha = p/(c_EE*N_E);
gamma = (1-rho_p)^4*0.028;  %integral calculated in mathematica
fun =@(x) 1/sqrt(2*pi)*exp(-x.^2/2).*max(A_EE*sqrt(alpha*gamma).*x+thrd,0);
omega_EE = integral(fun,-Inf,Inf);
g_I = sqrt(K_EE)*omega_EE/(J_EI*J_IE*sqrt(K_EI));
%% Creating EI/IE weight matrix
Jm_EI = (J_EI*cm_EI)/sqrt(K_EI);
Jm_IE = (J_IE*cm_IE)/K_IE;
%% Simulating network
% Familiar stimulus
I     = xi(:,4);
%I=rand(N_E,1);
r_E    = zeros(N_E,1)+sig(0);  %initial conditions
r_I    = zeros(N_I,1)+sig(0);  %initial conditions
flg = zeros(p-Pc,1);
ovlp  = zeros(Tlen,p-Pc);
for i = 1: Tlen
    if i<=t_stim_on          
        dr_E = (-r_E+J_EE*sig(r_E)-Jm_EI*g_I*max(r_I,0))/t_E*dt;
        dr_I = (-r_I+Jm_IE*sig(r_E))/t_I*dt;
    elseif i>t_stim_on && i <=t_stim_off
        dr_E = (-r_E+J_EE*sig(r_E+I)-Jm_EI*g_I*max(r_I,0))/t_E*dt;
        dr_I = (-r_I+Jm_IE*sig(r_E))/t_I*dt;
    else
        dr_E = (-r_E+J_EE*sig(r_E)-Jm_EI*g_I*max(r_I,0))/t_E*dt;
        dr_I = (-r_I+Jm_IE*sig(r_E))/t_I*dt;
    end
    r_E  = r_E+dr_E;
    r_I  = r_I+dr_I;
    rsamE(:,i) = r_E([rneuE]);
    rsamI(:,i) = r_I([rneuI]);
    for j=1:(p-Pc)
        ovlp(i,j) = corr(r_E,xi(:,j));
    end
end
%% Calculating correlations
nret = 0;  %number of retrieved patterns
corr_pat= zeros(p-Pc,1);
for i = 1:(p-Pc)
    corr_pat(i) =corr(r_E,sig(xi(:,i)));
    if corr_pat(i)>0.1
        nret = nret+1;
    end
end
%% Plotting Figures
close all
t_start = find(time==100);  %starting time is 100ms
x  = linspace(1,p-Pc,p-Pc);
figure
subplot(3,1,1)

for i =1:Pc
    plot(time(t_start:end)/1000,ovlp((t_start:end),(i-1)*(Pi-1)+1:i*(Pi-1),:),'linewidth',2)
    hold on
end
set(gca,'TickDir','out'); set(gca,'FontSize',14);
% text(0.2,0.9,str1,'Fontsize',12);
% text(0.5,0.9,str2,'Fontsize',12);
% text(0.75,0.9,str3,'Fontsize',12);
xline(t_stim_on/Tlen*tspan/1000); xline(t_stim_off/Tlen*tspan/1000);
xlabel('time (s)'); ylabel('Correlation'); 
title(['Correlation with familiar stimulus versus time when rho = ',num2str(rho)]);
ylim([0 1])

subplot(3,1,2)
for i =1:nsam
    plot(time(t_start:end)/1000,rsamE(i,t_start:end),'linewidth',2)
    hold on
end
ylim([-10 10])
set(gca,'TickDir','out'); set(gca,'layer','bottom');  set(gca,'FontSize',14);
xline(t_stim_on/Tlen*tspan/1000); xline(t_stim_off/Tlen*tspan/1000);
xlabel('time (s)'); ylabel('firing rate(Hz)');  
title('Firing rate of some excitatory neurons');
colorspec ={'red', 'blue','green','magenta','yellow'};
subplot(3,1,3)
for i =1:nsam
    plot(time(t_start:end)/1000,rsamI(i,t_start:end),'linewidth',2)
    hold on
end
ylim([0 1])
set(gca,'TickDir','out'); set(gca,'FontSize',18); set(gca,'FontSize',14);
xline(t_stim_on/Tlen*tspan/1000); xline(t_stim_off/Tlen*tspan/1000);
xlabel('time (s)'); ylabel('firing rate(Hz)');
title('Firing rate of some inhibitory neurons');
colorspec ={'red', 'blue','green','magenta','yellow'};
figure 
plot(x,corr_pat,'linewidth',2)
xlabel('memory'); ylabel('correlation'); axis tight
title('Correlation of final firing rate with different memories');
set(gca,'TickDir','out'); set(gca,'layer','bottom');  set(gca,'FontSize',14);

%%
save('Rate.mat','time','ovlp','rsamE','rsamI')
%% Miscellaneous Functions 
function y = sig(x)
    rm = 76.2;
    bt = 0.82;
    h0 = 2.46;
    y = rm./(1+exp(-bt*(x-h0)));
end
function y1 = f(x)
    xf = 26.6;
    bf = 0.28;
    qf = 0.83;
    y1 = 0.5*(2*qf-1 + tanh(bf*(x-xf)));
end
function y2 = g(x)
    xf = 26.6;
    bf = 0.28;
    qf = 0.99;
    y2 = 0.5*(2*qf-1 + tanh(bf*(x-xf)));
end