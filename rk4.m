function x1 = rk4(x0,h,f,par)
    
    k1 = h*f(x0,par);
    k2 = h*f(x0+k1/2,par);
    k3 = h*f(x0+k2/2,par);
    k4 = h*f(x0+k3,par);
    step = k1/6 + k2/3 + k3/3 + k4/6;
    Z = 10*sum((step.^2));
    x1 = x0 + step;%/Z;
end