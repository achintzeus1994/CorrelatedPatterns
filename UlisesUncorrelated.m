%% Parameter initialization
function ovlp_mem = UlisesUncorrelated(p,N)
tau  = 20;               % in ms 
tspan= 1000;             % simulation time= 2.5 sec
K     = 50000*0.005;      % number of synapses to each neuron
dt   = 0.5;              % time step 0.5 ms
A    = 3.55;
time = dt*linspace(0,tspan/dt,tspan/dt+1);
Tlen = length(time);
t_stim_on  = 0.4*Tlen;
t_stim_off = 0.6*Tlen;
net_dyn      = zeros(N,ceil(t_stim_off-t_stim_on));
%% Sample neurons
nsam = 20;
rneu = randi(N,nsam,1);
rsam = zeros(nsam,Tlen);
%% Recurrent matrix
nsyp = zeros(N,K);
for i =1:N
    nsyp(i,:) = randperm(N,K);    % number of input synapses/neurons
end
nneu = repmat(linspace(1,N,N)',1,K);
cm   = sparse(nneu,nsyp,1,N,N);
%% Creating patterns
% cormat = eye(p);
% L  = chol(cormat);  
% xi = null([normrnd(0,1,[N-(p+1),N]);ones(1,N)])*L;
% xi = xi - mean(xi);
% xi = xi./std(xi);
xi    = normrnd(0,1,[N,p]);
%% Creating weight matrix
xidiff_f = f(sig(xi));
xidiff_g = g(sig(xi));
J = A/K*cm.*(xidiff_f*xidiff_g');
for i = 1:N
    J(i,i)=0;  % Getting rid of autapses
end

%% Simulating network
% Familiar stimulus
mem_ret = 1;
I     = xi(:,mem_ret);
r    = zeros(N,1)+sig(0);  %initial conditions
ovlp  = zeros(Tlen,p);
for i = 1: Tlen
    if i<=t_stim_on       
        dr = (-r+sig(J*r))/tau*dt;       
    elseif i>t_stim_on && i <=t_stim_off
        net_dyn(:,i-ceil(t_stim_on)+1)=r;
        dr = (-r+sig(J*r+I))/tau*dt;
    else
        dr = (-r+sig(J*r))/tau*dt;
    end
    r  = r+dr;
    rsam(:,i) = r([rneu]);
    for j=1:p
       ovlp(i,j) = corr(r,g(sig(xi(:,j))));
    end
end
%% Calculating effective dimensionality
net_dyn=net_dyn(:,any(net_dyn));
C = cov(net_dyn);
e = eig(C);
lambda = e/sum(e);
dim   = 1/sum(lambda.^2);
%% Calculating correlations
nret = 0;  %number of retrieved patterns
corr_pat= zeros(p,1);
for i = 1:p
    corr_pat(i) =corr(r,g(sig(xi(:,i))));
    if corr_pat(i)>0.1
        nret = nret+1;
    end
end
ovlp_mem = corr_pat(mem_ret);
%% Plotting Figures
close all
t_start = find(time==100);  %starting time is 100ms
x  = linspace(1,p,p);
figure
subplot(2,1,1)
for i =1:min(p,5)
    plot(time(t_start:end)/1000,ovlp((t_start:end),:),'linewidth',2)
    hold on
end
set(gca,'TickDir','out'); set(gca,'layer','bottom');  set(gca,'FontSize',14);
str1 = 'Baseline'; str2 = 'Pres.'; str3 = 'Delay';
xline(t_stim_on/Tlen*tspan/1000); xline(t_stim_off/Tlen*tspan/1000);
text(0.2,0.9,str1,'Fontsize',12);
text(0.5,0.9,str2,'Fontsize',12);
text(0.75,0.9,str3,'Fontsize',12);
ylim([0 1])
xlabel('time (s)'); ylabel('Correlation');
subplot(2,1,2)
for i =1:nsam
    plot(time/1000,rsam(i,:),'linewidth',2)
    hold on
end
set(gca,'TickDir','out'); set(gca,'layer','bottom');  set(gca,'FontSize',14);
xline(t_stim_on/Tlen*tspan/1000); xline(t_stim_off/Tlen*tspan/1000);
xlabel('time (s)'); ylabel('firing rate(Hz)');
title('Firing rate of some neurons');
colorspec ={'red', 'blue','green','magenta','yellow'};
figure 
plot(x,corr_pat,'linewidth',2)
xlabel('memory'); ylabel('correlation'); axis tight
title('Correlation of final firing rate with different memories');

keyboard
end
