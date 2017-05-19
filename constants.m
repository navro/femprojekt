%% Load constants
qel = 9000;
alpha = 40;
Tinf = 20;
T0 = 30;

kvec = [66.8, 0.29, 1.059];
rhoc = [7265*210, 1850*950, 1850*950];
qnvec = [0, 0, Tinf*alpha, 0, 0, 0, qel, Tinf*alpha, 0, 0, 0, 0, 0];
nyvec = [0.36, 0.118, 0.136];
Evec = 1E9.*[50, 105, 105];
alphavec = [1.2E-5, 1.2E-5, 2E-5];