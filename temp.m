clear;
load halfmesh;
convert;
constants;
init_t;

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

% Plot the stationary solution
figure(1)
hold on
ed = extract(edof,astat);
fill(Ex', Ey', ed');
fill(-Ex', Ey', ed');

%% Calculate the C matrix
for i = 1:nelm
    Ce = plantml(Ex(i,:), Ey(i,:), rhoc(t(4,i)));
    C = assem(edof(i,:), C, Ce);
end

%% Time step
dt = 5;
cont = 1;

figure(2)
while 1 == cont
    a = solveq(K + Kc + C/dt, fb + C*a/dt);
    ed = extract(edof,a);
    clf
    hold on
    fill(Ex',Ey',ed');
    fill(-Ex',Ey',ed');
    colorbar;
    cont = 0;
    %waitforbuttonpress;
end

stress;