clear all; close all; clc;
%% Parameters
rho     = 0.25;
t_E     = 10;                       % time constant 
t_I     = 2;
N_E     = 10000;
N_I     = ceil(N_E/4);
N_tot   = N_E+N_I;
c_EE    = 2*0.025;         % EE connection probability
c_EI    = 0.02;
c_IE    = 0.02;
K_EE    = ceil(N_E*c_EE);
K_EI    = ceil(N_I*c_EI);
K_IE    = ceil(N_E*c_IE);
dt      = 0.5;                         % time bin in ms
tspan   = 1000;                       % total simulation time in ms
tauSE   = 20;
tauSI   = 5;
time    = dt*linspace(0,tspan/dt,tspan/dt+1);
Tlen    = length(time);
t_stim_on  = 0.3*Tlen;
t_stim_off = 0.7*Tlen;
Pc      = 5;        %number of classes
Pi      = 5;        %number of individual patterns/class
p       = Pc*Pi;
VrE  = -65;                   % reset voltage
VrI  = -65;
VthE = -50;                   % threshold voltage
VthI = -50;                   % threshold voltage
VminE = -75;                   % minimum voltage
VminI =-90;
V_E  = zeros(N_E,1)+VrE;     % Initial voltage of excitatory neurons in PCx
V_I  = zeros(N_I,1)+VrI;
t    = 0;
raster = zeros(N_tot,Tlen);
%% Calculating transfer function
close all;
VrE = 0;
VrI = 0;
s_noiseE = 2.7;
s_noiseI = 0.9;

VthE = 3.9;
VthI = 1.3;
%g_I = 278.77;
g_I = 304.595;
x   = linspace(-5,5,15);
nu_E= zeros(length(x),1);
nu_I= zeros(length(x),1);
fun = @(u) (1+erf(u)).*exp(u.^2);
sig_E = zeros(length(x),1);
sig_I = zeros(length(x),1);
for i=1: length(x)
    nu_E(i)  = 10^(-3)*t_E*sqrt(pi)*integral(fun,(VrE-x(i))/s_noiseE,(VthE-x(i))/s_noiseE);
    nu_I(i)  = 10^(-3)*t_I*sqrt(pi)*integral(fun,(VrI-x(i))/s_noiseI,(VthI-x(i))/s_noiseI);    
    sig_E(i) = sig(x(i));
    sig_I(i) = g_I*max(x(i),0);
end
nu_E=1./nu_E;
nu_I=1./nu_I;
figure
subplot(2,1,1)
plot(x,sig_E,x,nu_E,'linewidth',2);
legend('E Rate','Spike')
xlabel('x'); ylabel('nu');
title(sprintf('Excitatory rate/spiking transfer function ($V_r=$%d mV, $V_{th} =$%d mV, noise=%d mV)',VrE,VthE,s_noiseE),'Interpreter','latex')
axis tight
set(gca,'TickDir','out'); set(gca,'layer','bottom');  set(gca,'FontSize',14);
subplot(2,1,2)
plot(x,sig_I,x,nu_I,'linewidth',2);
legend('I Rate','Spike')
xlabel('x'); ylabel('nu');
title(sprintf('Inhibitory rate/spiking transfer function ($V_r=$%d mV, $V_{th} =$%d mV, noise=%d mV)',VrI,VthI,s_noiseI),'Interpreter','latex')
axis tight
set(gca,'TickDir','out'); set(gca,'layer','bottom');  set(gca,'FontSize',14);
deltaE = -65;
deltaI = -65;
VthstarE=-50;
VthstarI=-50;
lambdaE = (VthstarE-deltaE)/VthE;
lambdaI = (VthstarI-deltaI)/VthI;
%% New Parameters
VthE    = -50;
VthI    = -50;
VrE     = -65;
VrI     = -65;
s_noiseE= lambdaE*s_noiseE; 
s_noiseI= lambdaI*s_noiseI; 
%% Sample Neurons
nsam = 5;
VsamE = zeros(nsam,Tlen);
VsamI = zeros(nsam,Tlen);
%% Creating correlation matrix
cormat_cat = rho^2*(ones(Pi-1)-eye(Pi-1))+eye(Pi-1);   %desired correlation matrix
f_row  = rho*ones(1,Pi-1);
f_col  = [1; rho*ones(Pi-1,1)];
cormat_cat = [f_row;cormat_cat];
cormat_cat = [f_col, cormat_cat];
cormat     = cormat_cat;
for i =1:Pc-1
    cormat = blkdiag(cormat,cormat_cat);
end
% L  = chol(cormat);                         %Cholesky decomposition of correlation matrix
% xi = null([rand(N_E-(p+1),N_E);ones(1,N_E)])*L;
% xi = xi - mean(xi);
% xi = xi./std(xi);
% save('xi_rho_25_5pc_5pi.mat','xi'); 
xi=load('xi_rho_25_5pc_5pi.mat'); xi=xi.xi;
xic=xi(:,1);
for i=Pi:Pi:Pc*Pi-1
    xic=[xic,xi(:,i+1)];
end
del_var = 0:Pi:Pc*Pi-1;
xi(:,del_var+1)=[];
%% Creating weight matrix
[J_EE,Jm_EI,Jm_IE,g_I] = Weight_Matrix(N_E,N_I,K_EE,K_EI,K_IE,xi,xic,p,Pc,Pi,c_EE);
J_EE = lambdaE/tauSE*J_EE;
Jm_EI= lambdaE/tauSE*Jm_EI;
Jm_IE= lambdaI/tauSI*Jm_IE;
%% Network Dynamics
I    = sig(xi(:,2));
S_EE = zeros(N_E,1);
S_EI = zeros(N_E,1);
S_IE = zeros(N_I,1);
for i = 1:Tlen
    spk_EE = zeros(N_E,1);
    spk_EI = zeros(N_I,1);
    spk_IE = zeros(N_E,1);
    if i<=t_stim_on          
        V_E = V_E+dt/t_E*(VrE-V_E+S_EE-S_EI+s_noiseE*sqrt(t_E)*normrnd(0,1,[N_E,1]));
        V_I = V_I+dt/t_I*(VrI-V_I+S_IE+s_noiseI*sqrt(t_I)*normrnd(0,1,[N_I,1]));
    elseif i>t_stim_on && i <=t_stim_off
        V_E = V_E+dt/t_E*(VrE-V_E+S_EE-S_EI+I+s_noiseE*sqrt(t_E)*normrnd(0,1,[N_E,1]));
        V_I = V_I+dt/t_I*(VrI-V_I+S_IE+s_noiseI*sqrt(t_I)*normrnd(0,1,[N_I,1])); 
    else
        V_E = V_E+dt/t_E*(VrE-V_E+S_EE-S_EI+s_noiseE*sqrt(t_E)*normrnd(0,1,[N_E,1]));
        V_I = V_I+dt/t_I*(VrI-V_I+S_IE+s_noiseI*sqrt(t_I)*normrnd(0,1,[N_I,1]));
    end
    flagEth  = find(V_E>=VthE);
    flagEmin = find(V_E<VminE);
    flagIth  = find(V_I>=VthI);  
    flagImin = find(V_I<VminI);
    raster([flagEth;N_E+flagIth],i) = 1; 
    V_E(flagEth)  = VrE;
    V_E(flagEmin) = VminE;
    V_I(flagIth)  = VrI;
    V_I(flagImin) = VminI;
    spk_EE(flagEth) = 1;
    spk_EI(flagIth) = 1;
    spk_IE(flagEth) = 1;
    S_EE     = S_EE+dt/tauSE*(-S_EE+J_EE*tauSE*spk_EE);
    S_EI     = S_EI+dt/tauSE*(-S_EI+Jm_EI*tauSE*spk_EI);
    S_IE     = S_IE+dt/tauSI*(-S_IE+Jm_IE*tauSI*spk_IE);
    xint     = Jm_IE*tauSI*spk_IE;
    for j=1:nsam
        VsamE(j,i) = V_E(j);
        VsamI(j,i) = V_I(j);
    end
end
%%
for i =1:nsam
    flgE = find(VsamE(i,:)==-65);
    flgI = find(VsamI(i,:)==-65);
    VsamE(i,flgE)=10;
    VsamI(i,flgI)=10;
end


%% Rate calculation
fil_s=20;
x=-5*fil_s:dt:5*fil_s;
fun1 = @(u) 1/sqrt(2*pi)*exp(-u.^2/(2*fil_s^2));
x1   = fun1(x);
rate = conv2(raster,x1,'same');

r_E  = 1/(N_E)*sum(rate(1:N_E,:));                    %Calculating  E PSTH
r_I  = 1/(N_I)*sum(rate(N_E+1:N_tot,:));              %Calculating  E PSTH
a_E  = rate(1:N_E,:);
a_I  = rate(N_E+1:N_tot,:);
ovlp = zeros(p,Tlen);
for i=1:Tlen
    for j=1:Pc*(Pi-1)
        ovlp(j,i) = corr(a_E(:,i),sig(xi(:,j)));
    end
end
close all
figure
subplot(2,1,1)
for i =1:5
    plot(time/1000,a_E(i,:),'linewidth',2)
    hold on
end
xline(t_stim_on/Tlen*tspan/1000); xline(t_stim_off/Tlen*tspan/1000);
xlabel('time (s)'); ylabel('Excitatory Firing Rate(Hz)');
subplot(2,1,2)
for i =1:5
    plot(time/1000,a_I(i,:),'linewidth',2)
    hold on
end
xline(t_stim_on/Tlen*tspan/1000); xline(t_stim_off/Tlen*tspan/1000);
xlabel('time (s)'); ylabel('Inhibitory Firing Rate(Hz)');

%% Raster plot
% figure
% x=[0,time(end)]/1000;
% y=[0,1];
% imagesc(x,y,~raster(2,:))
% colormap('gray')
% xlabel('time (s)'); ylabel('Neurons');  
%xline(t_stim_on/Tlen*tspan/1000); xline(t_stim_off/Tlen*tspan/1000);

%% Plotting figures
figure
subplot(3,1,1)
for i = 1:p
    plot(time/1000,abs(ovlp(i,:)),'linewidth',2);
    hold on
end
xline(t_stim_on/Tlen*tspan/1000); xline(t_stim_off/Tlen*tspan/1000);
xlabel('time (s)'); ylabel('Correlation');
title(['Correlation with familiar stimulus versus time when rho = ',num2str(rho)]);
set(gca,'TickDir','out'); set(gca,'FontSize',14);
for i=1:nsam
    subplot(3,1,2)
    plot(time/1000,VsamE(i,:),'linewidth',2)
    xlabel('time (s)'); ylabel('Firing rate(Hz)');  
    title('Spiking of some excitatory neurons');
    xline(t_stim_on/Tlen*tspan/1000); xline(t_stim_off/Tlen*tspan/1000);
    set(gca,'TickDir','out'); set(gca,'FontSize',14);
    hold on
    subplot(3,1,3)
    plot(time/1000,VsamI(i,:),'linewidth',2)
    xlabel('time (s)'); ylabel('firing rate(Hz)');  
    title('Spiking of some inhibitory neurons');
    xline(t_stim_on/Tlen*tspan/1000); xline(t_stim_off/Tlen*tspan/1000);
    hold on
    set(gca,'TickDir','out'); set(gca,'FontSize',14);
end
figure
plot(time/1000,r_E,time/1000,r_I,'linewidth',2)
xline(t_stim_on/Tlen*tspan/1000); xline(t_stim_off/Tlen*tspan/1000);
xlabel('time (s)'); ylabel('Firing rate(Hz)');  
legend('Excitatory pop.','Inhibitory pop.')

%%
save('Spiking.mat','ovlp','VsamE','VsamI','time')