function xosav = Lorenz(xnp,par)

% sig = 10; 
% rho = 28; 
% beta = -8/3; 
% xnp = xnp(:);
x2 = [xnp(1,:).*xnp(2,:); xnp(1,:).*xnp(3,:); xnp(2,:).*xnp(3,:); xnp(1,:).^2; xnp(2,:).^2; xnp(3,:).^2];
x1 = [xnp; x2];
W = reshape(par,[3 9]);

%the Lorenz system

% %Test function 1
xosav = W*x1;

end
   