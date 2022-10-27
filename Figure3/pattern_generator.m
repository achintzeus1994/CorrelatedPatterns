function xi = pattern_generator(lp_no,rho)
p1 = lp_no-1;
rho1= int8(10*rho);
dat= load(sprintf('fp_32_rho_%d_N_10000_p_val_%d.mat',rho1,p1));
xi =dat.xi;
xi =double(xi);
end