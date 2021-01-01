clear all; close all; clc;
rng(2);
N     = 10000;
alpha = linspace(0.1,0.1,1);
K     = 50000*0.005;      % number of synapses to each neuron
p     = ceil(alpha*K);
c     = K/N;              % connection probability
len   = length(p);
ovlp_mem = zeros(len,1);
for i = 1:len
    ovlp_mem(i) = UlisesUncorrelated(p(i),N);
end

figure
plot(alpha,ovlp_mem,'linewidth',2)
%%
%corr_sim00 = [alpha',ovlp_mem];
%save('corr_sim00.mat','corr_sim00')