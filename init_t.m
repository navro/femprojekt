%% Initialize matrices for the temperature part of the project
K = zeros(ndof);
Kc = zeros(ndof);
C = zeros(ndof);
a = T0.*ones(ndof,1);
fb = zeros(ndof,1);