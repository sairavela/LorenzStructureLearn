
dt = 0.01; %time step

sig = 10;
rho = 28; 
beta = 8/3;

W = [-sig sig 0 0 0 0 0 0 0; rho -1 0 0 -1 0 0 0 0; 0 0 -beta 1 0 0 0 0 0]*dt; %Euler scheme
x0 = [-1.1; 2.2; -2.7]; %initial condition
par = [W(:)]; % parameters
