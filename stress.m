%% Initialize matrices
K = zeros(2*ndof);
f0 = zeros(2*ndof,1);
newedof = [edof(:,1), edof(:,2), edof(:,2) + ndof, edof(:,3), edof(:,3) + ndof, edof(:,4), edof(:,4) + ndof];
seff_el = zeros(nelm,1);
seff_nod = zeros(ndof,1);

%% Calculate the K matrix and f0
for i = 1:nelm
        E = Evec(t(4,i));
        ny = nyvec(t(4,i));
        D = E/((1+ny)*(1-2*ny)) * [1-ny, ny, 0; ny, 1-ny, 0; 0, 0, (1-2*ny)/2];
        dt = mean(astat(t(1:3, i)));
        e0 = (1 + ny) * alphavec(t(4,i)) * dt * [1; 1; 0];
        
        % K matrix
        Ke = plante(Ex(i,:), Ey(i,:), [2, 1], D);
        
        % f0 matrix
        f0e = plantf(Ex(i,:), Ey(i,:), [2, 1], (D*e0)');
        
        % Assemble!!
        [K,f0] = assem(newedof(i,:), K, Ke, f0, f0e);
end

%% Boundary conditions
bc = [1:2*ndof; zeros(1,2*ndof); zeros(1,2*ndof)]'; 
for i = 1:length(e)
    if e(5,i) == 1
        bc(e(1,i) + ndof,3) =  1;
        bc(e(2,i) + ndof,3) =  1;
    elseif e(5,i) == 2 || e(5,i) == 12
        bc(e(1,i),3) =  1;
        bc(e(2,i),3) =  1;
    end
 end
 
 % Remove rows without a BC
 empty = find(bc(:,3) == 0);
 bc(empty,:) = [];
 bc = bc(:,1:2);
 
 %% Solve the equation!
 [a,r] = solveq(K,f0,bc);
 
 %% Compute the von Mises stress for each element
 for i = 1:nelm
     ed = extract(newedof,a);
     E = Evec(t(4,i));
     ny = nyvec(t(4,i));
     D = E/((1+ny)*(1-2*ny)) * [1-ny, ny, 0; ny, 1-ny, 0; 0, 0, (1-2*ny)/2];
    [es,et] = plants(Ex(i,:), Ey(i,:), [2, 1], D, ed);
     
    tau = E/(2*(1+ny)) * et(4:6);
    seff_el(i) = sqrt(es(1)^2 + es(2)^2 + es(3)^2 + es(1)*es(2) + es(1)*es(3) + es(2)*es(3) + 3*tau(1)^2 + 3*tau(2)^2 + 3*tau(3)^2); 
 end
 
 %% Compute the von Mises stress for each node
for i=1:size(coord,1)
    [c0,c1]=find(edof(:,2:4)==i);
    seff_nod(i,1)=sum(seff_el(c0))/size(c0,1);
end

%% Plot the stresses
hold on
fill(Ex', Ey', seff_el');
fill(-Ex', Ey', seff_el');