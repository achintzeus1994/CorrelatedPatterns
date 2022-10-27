function ovlp_mem = Gutfreund_Rate(lp_no,rho,Pc,Pi,xi_ret,N,c)
%% Parameter initialization
tau  = 20;               % in ms 
K    = N*c;              % number of synapses to each neuron
tspan= 1000;             % simulation time= 2.5 sec
dt   = 0.5;              % time step 0.5 ms
A    = 4;
time = dt*linspace(0,tspan/dt,tspan/dt+1);
Tlen = length(time);
t_stim_on  = 0.4*Tlen;
t_stim_off = 0.6*Tlen;
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
%% Generating patterns
xi = pattern_generator(lp_no,rho);
xic=xi(:,1);
del_var = 0:Pi:Pc*Pi-1;
xi(:,del_var+1)=[];
Pi = size(xi,2);
p   = Pc*Pi;
%% Calculating pattern correlations
rho_f = zeros(Pi,1);
rho_g = zeros(Pi,1);
flg   = 1;
for j =1:Pi
    rho_f(flg) = corr(f(sig(xic)),f(sig(xi(:,j))));
    rho_g(flg) = corr(g(sig(xic)),g(sig(xi(:,j))));
    flg = flg + 1;
end
rho_f = mean(rho_f);
rho_g = mean(rho_g);
%% Creating weight matrix
xidiff_f = f(sig(xi))-rho_f*f(sig(repelem(xic,1,Pi)));
xidiff_g = g(sig(xi))-rho_g*g(sig(repelem(xic,1,Pi)));
J = A/K*cm.*(xidiff_f*xidiff_g');
for i = 1:N
    J(i,i)=0;  % Getting rid of autapses
end
%% Simulating network
% Familiar stimulus
mem_ret=xi_ret;
I     = xi(:,mem_ret);
r    = zeros(N,1)+sig(0);  %initial conditions
ovlp  = zeros(Tlen,Pi);
for i = 1: Tlen
    if i<=t_stim_on       
        dr = (-r+sig(J*r))/tau*dt;       
    elseif i>t_stim_on && i <=t_stim_off
        dr = (-r+sig(J*r+I))/tau*dt;
    else
        dr = (-r+sig(J*r))/tau*dt;
    end
    r  = r+dr;
    rsam(:,i) = r([rneu]);
    for j=1:Pi
       ovlp(i,j) = corr(r,g(sig(xi(:,j))));
    end
end
%% Calculating correlations
nret = 0;  %number of retrieved patterns
corr_pat= zeros(Pi,1);
for i = 1:(p-Pc)
    corr_pat(i) =corr(r,g(sig(xi(:,i))));
    if corr_pat(i)>0.1
        nret = nret+1;
    end
end
ovlp_mem = corr_pat(mem_ret);
%% Plotting Figures
close all
x  = linspace(1,p,p);
colorspec ={'red', 'blue','green','magenta','yellow','cyan'};
figure

subplot(2,1,1)
for i = 1:Pi
    plot(time/1000,abs(ovlp(:,i)),'linewidth',2,'Color', colorspec{1})
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
title(['Correlation with familiar stimulus versus time when rho = ',num2str(rho)]);

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
end