clear;
load halfmesh;
convert;
constants;
initialize;

%% Half mesh values
kvec = [66.8, 0.29, 1.059];
rhoc = [7265*210, 1850*950, 1850*950];
qnvec = [0, 0, Tinf*alpha, 0, 0, 0, qel, Tinf*alpha, 0, 0, 0, 0, 0];

%% Calculate the K matrix
for i = 1:nelm
    k = kvec(t(4,i));
    Ke = flw2te(Ex(i,:), Ey(i,:), 1, k);
    K = assem(edof(i,:), K, Ke);
end

%% Calculate the boundary vector and the Kc matrix
for i = 1:length(e)
    seg = e(5,i);
    n1 = e(1,i);
    n2 = e(2,i);
    L = sqrt((p(1,n1) - p(1,n2))^2 + (p(2,n1) - p(2,n2))^2);
    
    fb(n1) = fb(n1) + L/2*qnvec(seg);
    fb(n2) = fb(n2) + L/2*qnvec(seg);
    
    if seg == 3 || seg == 8
        Kce = L*alpha/6 * [2, 1; 1, 2];
        Kc = assem([0, n1, n2], Kc, Kce);
    end
end

%% Stationary solution
astat = solveq(K + Kc, fb);

%% Calculate the C matrix
for i = 1:nelm
    Ce = plantml(Ex(i,:), Ey(i,:), rhoc(t(4,i)));
    C = assem(edof(i,:), C, Ce);
end

%% Time step
dt = 5;

while 1 == 1
    a = solveq(K + Kc + C/dt, fb + C*a/dt);
    ed = extract(edof,a);
    clf
    hold on
    fill(Ex',Ey',ed');
    fill(-Ex',Ey',ed');
    colorbar;
    waitforbuttonpress;
end