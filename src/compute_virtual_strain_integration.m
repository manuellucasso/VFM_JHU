function [delta_e, grad_delta_u] = compute_virtual_strain_integration(model, coords, delta_u, gauss_order)
% compute_virtual_strain
% Computes virtual strain tensor delta(e) = sym(grad(delta(u)))
% for all elements at specified Gauss points.
%
% INPUTS:
%   model       : struct with .elements (nElem x 9) connectivity
%   coords      : (nNodes x 3) nodal coordinates
%   delta_u     : (nNodes x 3) displacement difference (u_guess - u_truth)
%   gauss_order : integration order (1, 2, 3, etc.)
%
% OUTPUT:
%   delta_e       : (nElem x nGauss x 3 x 3) virtual strain tensor at each Gauss point
%   grad_delta_u  : (nElem x nGauss x 3 x 3) displacement gradient at each Gauss point
%   gauss_points  : (nGauss x 3) Gauss point coordinates [xi, eta, zeta]

nElem = size(model.elements, 1);

% Get Gauss points for hexahedral element
[gauss_points, ~] = get_gauss_points(3, gauss_order);
nGauss = size(gauss_points, 1);

% Preallocate outputs
delta_e = zeros(nElem, nGauss, 3, 3);
grad_delta_u = zeros(nElem, nGauss, 3, 3);

for k = 1:nElem
    % Get nodal indices for this element
    nodes = model.elements(k, 2:9);
    
    % Get nodal coordinates
    X = coords(nodes, :);  % 8x3 matrix
    
    % Get virtual displacements for these nodes
    delta_u_nodes = delta_u(nodes, :);  % 8x3 matrix
    
    for gp = 1:nGauss
        xi = gauss_points(gp, 1);
        eta = gauss_points(gp, 2);
        zeta = gauss_points(gp, 3);
        
        % Calculate shape function derivatives at this Gauss point
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
        
        % Reference tangents
        dX_dxi = dN_dxi' * X;    % 1x3 vector
        dX_deta = dN_deta' * X;  % 1x3 vector
        dX_dzeta = dN_dzeta' * X; % 1x3 vector
        
        % Jacobian matrix
        jacmat = [dX_dxi', dX_deta', dX_dzeta'];  % 3x3
        
        % Shape function derivatives in physical coordinates
        dN_dX = [dN_dxi, dN_deta, dN_dzeta] / jacmat;  % 8x3 matrix
        
        % Calculate grad(delta_u) = sum_i delta_u_i ⊗ ∇N_i
        gdu = zeros(3, 3);
        for i = 1:8
            gdu = gdu + delta_u_nodes(i, :)' * dN_dX(i, :);
        end
        
        % Save results
        grad_delta_u(k, gp, :, :) = gdu;
        delta_e(k, gp, :, :) = 0.5 * (gdu + gdu');
    end
end
end