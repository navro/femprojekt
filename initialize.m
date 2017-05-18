%% Initialize matrices
K = zeros(ndof);
Kc = zeros(ndof);
Kc1 = zeros(ndof);
C = zeros(ndof);
a = T0.*ones(ndof,1);
fb = zeros(ndof,1);