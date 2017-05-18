clear;
load mesh;
convert;
constants;
initialize;

%% Calculate the K matrix
for i = 1:nelm
    k = kvec(t(4,i));
    Ke = flw2te(Ex(i,:), Ey(i,:), 1, k);
    K = assem(edof(i,:), K, Ke);
end

%% Calculate the boundary vector
for i = 1:length(e)
    seg = e(5,i);
    res = 0;
    
    if seg == 2 || seg == 4 || seg == 13 || seg == 15   % Convection
        res = Tinf*alpha;
    elseif seg == 14                                    % "Heat Source"
        res = qel;
    else
        continue
    end
    
    n1 = e(1,i); n2 = e(2,i);
    L = sqrt((p(1,n1) - p(1,n2))^2 + (p(2,n1) - p(2,n2))^2);
    
    fb(n1) = fb(n1) + L*res/2;
    fb(n2) = fb(n2) + L*res/2;
end

%% Calculate the Kc matrix (from the convection)
for i = 1:length(e)
    seg = e(5,i);
    
    if seg == 4 || seg == 13 || seg == 2 || seg == 15
        n1 = e(1,i);
        n2 = e(2,i);
        dx = p(1, n1) - p(1,n2);
        dy = p(2, n1) - p(2,n2);
        L = sqrt(dx.^2 + dy.^2);
        
        Kce = L*alpha/6* [2, 1; 1, 2];
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
dt = 1;

while 1 == 1
    a = solveq(K + Kc + C/dt, fb + C*a/dt);
    ed = extract(edof,a);
    clf
    fill(Ex',Ey',ed');
    colorbar;
    waitforbuttonpress;
end