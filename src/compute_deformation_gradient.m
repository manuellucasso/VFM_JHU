function [F_all, J_all] = compute_deformation_gradient(model, coords, disp, gauss_order)
% Compute deformation gradient F at each Gauss point
% model.elements: nElem x 9 connectivity (cols 2:9 = 8 node indices)
% coords: (nNodes x 3) reference coordinates
% disp:   (nNodes x 3) displacements
% gauss_order: integration order

nElem = size(model.elements,1);
[gauss_points, weights] = get_gauss_points(3, gauss_order);
nGauss = size(gauss_points,1);

% Allocate
F_all = zeros(nElem, nGauss, 3, 3);
J_all = zeros(nElem, nGauss);

% Current configuration
xcoords = coords + disp;

for e = 1:nElem
    nodes = model.elements(e,2:9);  % 8-node hex
    X = coords(nodes,:);   % reference
    x = xcoords(nodes,:);  % deformed
    
    for gp = 1:nGauss
        xi   = gauss_points(gp,1);
        eta  = gauss_points(gp,2);
        zeta = gauss_points(gp,3);
        
        % --- shape function derivatives wrt xi,eta,zeta ---
        dN_dxi = 0.125 * [-(1-eta)*(1-zeta);  (1-eta)*(1-zeta); ...
                           (1+eta)*(1-zeta); -(1+eta)*(1-zeta); ...
                          -(1-eta)*(1+zeta);  (1-eta)*(1+zeta); ...
                           (1+eta)*(1+zeta); -(1+eta)*(1+zeta)];
        
        dN_deta = 0.125 * [-(1-xi)*(1-zeta); -(1+xi)*(1-zeta); ...
                            (1+xi)*(1-zeta);  (1-xi)*(1-zeta); ...
                           -(1-xi)*(1+zeta); -(1+xi)*(1+zeta); ...
                            (1+xi)*(1+zeta);  (1-xi)*(1+zeta)];
        
        dN_dzeta = 0.125 * [-(1-xi)*(1-eta); -(1+xi)*(1-eta); ...
                            -(1+xi)*(1+eta); -(1-xi)*(1+eta); ...
                             (1-xi)*(1-eta);  (1+xi)*(1-eta); ...
                             (1+xi)*(1+eta);  (1-xi)*(1+eta)];
        
        % Jacobian wrt reference element
        dX_dxi   = dN_dxi'   * X;
        dX_deta  = dN_deta'  * X;
        dX_dzeta = dN_dzeta' * X;
        
        Jmat = [dX_dxi', dX_deta', dX_dzeta'];   % 3x3 Jacobian
        dN_dX = [dN_dxi dN_deta dN_dzeta] / Jmat;  % (8x3)
        
        % --- Compute F ---
        F = zeros(3,3);
        for i=1:8
            F = F + x(i,:)' * dN_dX(i,:);   % outer product
        end
        
        F_all(e,gp,:,:) = F;
        J_all(e,gp) = det(F);
    end
end
end
 