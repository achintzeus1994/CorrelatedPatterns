function y = sig(x)
    rm = 76.2;
    bt = 0.82;
    h0 = 2.46;
    y = rm./(1+exp(-bt*(x-h0)));
end
