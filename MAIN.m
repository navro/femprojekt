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
alphaSolder = 1.2E-5;
alphaSMD = 1.2E-5;
Tinf = 20;
T0 = 30;

qn = [0, Tinf*alphaSolder, 0, Tinf*alphaSolder, -37, -37, 0, 0, -37, -37, -37, -37, Tinf*alphaSMD, qel,  Tinf*alphaSMD, -37, -37, -37, -37, -37, -37]; 

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
