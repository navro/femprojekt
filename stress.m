init_s;

%% Create a new edof
% We need to create a new edof since we have two values (ux, uy) at each
% node.
newedof = [edof(:,1), edof(:,2), edof(:,2) + ndof, edof(:,3), edof(:,3) + ndof, edof(:,4), edof(:,4) + ndof];

%% Calculate the K matrix and f0
for i = 1:nelm
        E = Evec(t(4,i));   % Defined in constanta
        ny = nyvec(t(4,i)); % Defined in constants
        D = E/((1+ny)*(1-2*ny)) * [1-ny, ny, 0; ny, 1-ny, 0; 0, 0, (1-2*ny)/2];
        
        % Approximate the element's temperature with the mean value of the
        % temperature at the nodes
        dt = mean(astat(t(1:3, i))) - T0; 
        e0 = (1 + ny) * alphavec(t(4,i)) * dt * [1; 1; 0];
        
        % Generate the element matrices Ke and f0e
        Ke = plante(Ex(i,:), Ey(i,:), [2, 1], D);
        f0e = plantf(Ex(i,:), Ey(i,:), [2, 1], (D*e0)');
        
        % Assemble!!
        [K,f0] = assem(newedof(i,:), K, Ke, f0, f0e);
end

%% Boundary conditions
% Let the first column represent the dof index, and the second the node's
% boundary condition. The third column is temporary and marks if the row
% should be included in the final matrix.
bc = [1:2*ndof; zeros(1,2*ndof); zeros(1,2*ndof)]';

% Loop through the edge lines
for i = 1:length(e)
    if e(5,i) == 1  % uy = 0
        bc(e(1,i) + ndof,3) = 1;    % Mark edge node 1 for inclusion (y component only)
        bc(e(2,i) + ndof,3) = 1;    % Mark edge node 2 for inclusion (y component only)
    elseif e(5,i) == 2 || e(5,i) == 12  || e(5,i) == 13  % ux = 0
        bc(e(1,i),3) =  1;          % Mark edge node 1 for inclusion (x component only)
        bc(e(2,i),3) =  1;          % Mark edge node 2 for inclusion (x component only)
    end
 end
 
 % Remove the rows that do not have a BC (where the last column is 0)
 empty = find(bc(:,3) == 0);
 bc(empty,:) = [];
 
 % Delete the last colum to get the format that CALFEM expects
 bc = bc(:,1:2);
 
 %% Solve the equation!
 [a,r] = solveq(K,f0,bc);
 
 %% Compute the von Mises stress for each element
 for i = 1:nelm
    % Extract displacements for the element
    ed = extract(newedof(i,:),a);
    
    % Load constants & temperature approximation
    E = Evec(t(4,i));
    ny = nyvec(t(4,i));
    D = E/((1+ny)*(1-2*ny)) * [1-ny, ny, 0; ny, 1-ny, 0; 0, 0, (1-2*ny)/2];
    dt = mean(astat(t(1:3, i))) - T0;
    
    % Let CALFEM calculate the stresses and strains for the element
    [es,et] = plants(Ex(i,:), Ey(i,:), [2, 1], D, ed);
    
    % Calculate sigma_zz manually (not available in es)
    sig_zz = ny*(es(1)+es(2)) - alphavec(t(4,i)) * E * dt;
    
    % Calculate tau_xy from gamma_xy
    tau_xy = E/(2*(1+ny)) * et(3);
    
    % Subtract thermoelastic strains from es
    es = es - D*e0;
    
    % Calculate the von Mises stress for the element
    seff_el(i) = sqrt(es(1)^2 + es(2)^2 + sig_zz^2 - es(1)*es(2) - es(1)*sig_zz - es(2)*es(3) + 3*tau_xy); 
 end
 
 %% Compute the von Mises stress for each node
for i=1:size(coord,1)
    [c0,c1]=find(edof(:,2:4)==i);
    seff_nod(i,1)=sum(seff_el(c0))/size(c0,1);
end

%% Plot the stresses
figure(3)
hold on
seff_eff = extract(edof, seff_nod);
fill(Ex', Ey', seff_eff');
fill(-Ex', Ey', seff_eff');

%% Plot the displacements
figure(4)
hold on
us = extract(newedof, a);
sfac = eldisp2(Ex, Ey, us,[1 4 0],1000);