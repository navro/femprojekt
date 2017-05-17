%% Convert the mesh to the CALFEM format
nelm=length(t(1,:));
edof(:,1)=1:nelm;
edof(:,2:4)=t(1:3,:)';
coord=p';
ndof=max(max(t(1:3,:)));
[Ex,Ey]=coordxtr(edof,coord,(1:ndof)',3);

%% Initialize and assign values
kvec = [66.8, 66.8, 1.059, 0.29];
K = zeros(ndof);
D = eye(2);
ep = 1;   % Thickness
qel = 9000;
alpha = 1.2E-5;
alpha = 1.2E-5;
Tinf = 20;
T0 = 30;

qn = [0, Tinf*alpha, 0, Tinf*alpha, -37, -37, 0, 0, -37, -37, -37, -37, Tinf*alpha, qel,  Tinf*alpha, -37, -37, -37, -37, -37, -37]; 
rhoc = [7265*210, 7265*210, 1850*950, 1850*950];
d0 = T0.*ones(ndof,1);

%% Calculate the K matrix
for i = 1:nelm
    k = kvec(t(4,i));
    Ke = flw2te(Ex(i,:), Ey(i,:), ep, k.*D);
    %K = assem(edof, K, Ke);
    K = assem(edof(i,:), K, Ke);
end

%% Calculate the boundary vector
bc = [1:ndof; zeros(1,ndof); zeros(1,ndof)]'; 
for i = 1:length(e)
    if qn(e(5,i)) ~= -37
        bc(e(1,i),2) =  bc(e(1,i),2) + qn(e(5,i))./2;
        bc(e(2,i),2) =  bc(e(2,i),2) + qn(e(5,i))./2;
        bc(e(1,i),3) =  1;
        bc(e(2,i),3) =  1;
    end
end

% Remove internal rows
empty = find(bc(:,3) == 0);
bc(empty,:) = [];

%% Calculate the Kc matrix (from the convection)
Kc = zeros(ndof);
J = 1/6 .* [1 2; 2 1];

for i = 1:length(e)
    if e(5,i) == 4 || e(5,i) == 13 || e(5,i) == 2 || e(5,i) == 15
        dx = p(1, e(1,i)) - p(1, e(2,i));
        dy = p(2, e(1,i)) - p(2, e(2,i));
        L = sqrt(dx^2 + dy^2);
        Ne = L.*J;
        Kc(e(1,i), e(1,i)) = Kc(e(1,i), e(1,i)) + L.*alpha./6;
        Kc(e(2,i), e(2,i)) = Kc(e(2,i), e(2,i)) + L.*alpha./6;
        Kc(e(1,i), e(2,i)) = Kc(e(1,i), e(2,i)) + L.*alpha./3;
        Kc(e(2,i), e(1,i)) = Kc(e(2,i), e(1,i)) + L.*alpha./3;
    end
end

Ktot = K + Kc;
%a = solveq(Ktot,zeros(ndof,1),bc);

%% Calculate the Ce matrix
Kc = zeros(ndof);
for i = 1:nelm
    Kce = plantml(Ex(i,:), Ey(i,:), rhoc(t(4,i)));
    Kc(t(1,i),t(1,i)) = Kc(t(1,i),t(1,i)) + Kce(1,1);
    Kc(t(1,i),t(2,i)) = Kc(t(1,i),t(2,i)) + Kce(1,2);
    Kc(t(1,i),t(3,i)) = Kc(t(1,i),t(3,i)) + Kce(1,3);
    
    Kc(t(2,i),t(1,i)) = Kc(t(2,i),t(1,i)) + Kce(2,1);
    Kc(t(2,i),t(2,i)) = Kc(t(2,i),t(2,i)) + Kce(2,2);
    Kc(t(2,i),t(3,i)) = Kc(t(2,i),t(3,i)) + Kce(2,3);
    
    Kc(t(3,i),t(1,i)) = Kc(t(3,i),t(1,i)) + Kce(3,1);
    Kc(t(3,i),t(2,i)) = Kc(t(3,i),t(2,i)) + Kce(3,2);
    Kc(t(3,i),t(3,i)) = Kc(t(3,i),t(3,i)) + Kce(3,3);
end

%% Time dependent things
dt = 0.02;
T = 5;
alpha = 0.5;    % Trapezoidal rule = implicit method