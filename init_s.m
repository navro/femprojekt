%% Initialize matrices for the stress part of the project
K = zeros(2*ndof);
f0 = zeros(2*ndof,1);
seff_el = zeros(nelm,1);
seff_nod = zeros(ndof,1);