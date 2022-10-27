function y2 = g(x)
    xf = 26.6;
    bf = 0.28;
    qf = 0.99;
    y2 = 0.5*(2*qf-1 + tanh(bf*(x-xf)));
end