function y1 = f(x)
    xf = 26.6;
    bf = 0.28;
    qf = 0.83;
    y1 = 0.5*(2*qf-1 + tanh(bf*(x-xf)));
end

