clear all; close all;clc;
Pc  = 5;        %number of classes
Pi  = 5;        %number of individual patterns/class
p   = Pc*Pi;
rho = 0.25; 
%% Rate
rate_Un = load('Rate_Un.mat');
rsamUn  = rate_Un.rsam;
ovlpUn  = rate_Un.ovlp;

%% Rate EI
rateEI = load('Rate.mat');
EIrate_ovlp = rateEI.ovlp;
rsamE= rateEI.rsamE;
rsamI= rateEI.rsamI;
time = rateEI.time;
Tlen = length(time);
tspan= 1000;  
nsam = 20;
t_stim_on  = 0.4*Tlen;
t_stim_off = 0.6*Tlen;
%% Spiking EI
spikeEI = load('Spiking.mat');
EIspike_ovlp = (spikeEI.ovlp)';
VsamE= spikeEI.VsamE;
VsamI= spikeEI.VsamI;
%%
close all
t_start = find(time==100);  %starting time is 100ms
x  = linspace(1,p-Pc,p-Pc);
figure
%%
subplot(3,3,1)
for i =1:Pc
    plot(time(t_start:end)/1000,ovlpUn((t_start:end),(i-1)*(Pi-1)+1:i*(Pi-1),:),'linewidth',2)
    hold on
end
axis tight
set(gca,'TickDir','out'); set(gca,'FontSize',14);
xline(t_stim_on/Tlen*tspan/1000); xline(t_stim_off/Tlen*tspan/1000);
xlabel('time (s)'); ylabel('Correlation'); 
title([{'Correlation of familiar  '},{'stimulus versus time rho = 0.25'}]);
ylim([0 1])
%%
subplot(3,3,2)

for i =1:Pc
    plot(time(t_start:end)/1000,EIrate_ovlp((t_start:end),(i-1)*(Pi-1)+1:i*(Pi-1),:),'linewidth',2)
    hold on
end
axis tight
set(gca,'TickDir','out'); set(gca,'FontSize',14);
xline(t_stim_on/Tlen*tspan/1000); xline(t_stim_off/Tlen*tspan/1000);
xlabel('time (s)'); ylabel('Correlation'); 
title([{'Correlation of familiar  '},{'stimulus versus time rho = 0.25'}]);
ylim([0 1])
%%
subplot(3,3,3)

for i =1:Pc
    plot(time(t_start:end)/1000,EIspike_ovlp((t_start:end),(i-1)*(Pi-1)+1:i*(Pi-1),:),'linewidth',2)
    hold on
end
axis tight
set(gca,'TickDir','out'); set(gca,'FontSize',14);
xline(t_stim_on/Tlen*tspan/1000); xline(t_stim_off/Tlen*tspan/1000);
xlabel('time (s)'); ylabel('Correlation'); 
title([{'Correlation of familiar  '},{'stimulus versus time rho = 0.25'}]);
ylim([0 1])
%%
subplot(3,3,4)

for i =1:nsam
    plot(time(t_start:end)/1000,rsamUn(i,t_start:end),'linewidth',2)
    hold on
end
axis tight
set(gca,'TickDir','out'); set(gca,'layer','bottom');  set(gca,'FontSize',14);
xline(t_stim_on/Tlen*tspan/1000); xline(t_stim_off/Tlen*tspan/1000);
xlabel('time (s)'); ylabel('firing rate(Hz)');  
title([{'Firing rate of'},{'some neurons'}]);
colorspec ={'red', 'blue','green','magenta','yellow'};

%%
subplot(3,3,5)
for i =1:nsam
    plot(time(t_start:end)/1000,rsamE(i,t_start:end),'linewidth',2)
    hold on
end
axis tight
ylim([-10 10])
set(gca,'TickDir','out'); set(gca,'layer','bottom');  set(gca,'FontSize',14);
xline(t_stim_on/Tlen*tspan/1000); xline(t_stim_off/Tlen*tspan/1000);
xlabel('time (s)'); ylabel('firing rate(Hz)');  
title([{'Firing rate of'},{'some excitatory neurons'}]);
colorspec ={'red', 'blue','green','magenta','yellow'};
%%
subplot(3,3,8)
for i =1:nsam
    plot(time(t_start:end)/1000,rsamI(i,t_start:end),'linewidth',2)
    hold on
end
axis tight
ylim([0 1])
set(gca,'TickDir','out'); set(gca,'FontSize',18); set(gca,'FontSize',14);
xline(t_stim_on/Tlen*tspan/1000); xline(t_stim_off/Tlen*tspan/1000);
xlabel('time (s)'); ylabel('firing rate(Hz)');
title([{'Firing rate of'},{'some inhibitory neurons'}]);
colorspec ={'red', 'blue','green','magenta','yellow'};
%%
for i=1:nsam
    subplot(3,3,6)
    plot(time/1000,VsamE(i,:),'linewidth',2)
    xlabel('time (s)'); ylabel('Firing rate(Hz)');  
    title([{'Spiking of'},{'some excitatory neurons'}]);
    xline(t_stim_on/Tlen*tspan/1000); xline(t_stim_off/Tlen*tspan/1000);
    set(gca,'TickDir','out'); set(gca,'FontSize',14);
    hold on
    subplot(3,3,9)
    plot(time/1000,VsamI(i,:),'linewidth',2)
    xlabel('time (s)'); ylabel('firing rate(Hz)');  
    title([{'Spiking of'},{'some inhibitory neurons'}]);
    xline(t_stim_on/Tlen*tspan/1000); xline(t_stim_off/Tlen*tspan/1000);
    hold on
    set(gca,'TickDir','out'); set(gca,'FontSize',14);
end