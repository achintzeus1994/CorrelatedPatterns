function [J_EE,Jm_EI,Jm_IE,g_I] =Weight_Matrix(N_E,N_I,K_EE,K_EI,K_IE,xi,xic,p,Pc,Pi,c_EE)
%% Recurrent matrix
A_EE = 1.5*4;
rho_p = zeros(p-Pc,1);
flg   = 1;
for i=1:Pc
    for j =1:Pi-1
        rho_p(flg) = corr(f(sig(xic(:,i))),f(sig(xi(:,j+(i-1)*(Pi-1)))));
        flg = flg + 1;
    end
end
rho_p = mean(rho_p);
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
end