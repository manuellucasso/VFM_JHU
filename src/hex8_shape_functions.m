function [N, dN_dxi, dN_deta, dN_dzeta] = hex8_shape_functions(xi_eta_zeta)
%HEX8_SHAPE_FUNCTIONS Evaluates the 8-node hex shape functions and their derivatives
% Input:
%   xi_eta_zeta: [xi, eta, zeta] (1x3)
% Output:
%   N:        [8x1] shape function values
%   dN_dxi:   [8x1] derivatives wrt xi
%   dN_deta:  [8x1] derivatives wrt eta
%   dN_dzeta: [8x1] derivatives wrt zeta

xi   = xi_eta_zeta(1);
eta  = xi_eta_zeta(2);
zeta = xi_eta_zeta(3);

N = 1/8 * [
    (1-xi)*(1-eta)*(1-zeta);
    (1+xi)*(1-eta)*(1-zeta);
    (1+xi)*(1+eta)*(1-zeta);
    (1-xi)*(1+eta)*(1-zeta);
    (1-xi)*(1-eta)*(1+zeta);
    (1+xi)*(1-eta)*(1+zeta);
    (1+xi)*(1+eta)*(1+zeta);
    (1-xi)*(1+eta)*(1+zeta);
];

dN_dxi = 1/8 * [
    -(1-eta)*(1-zeta);
     (1-eta)*(1-zeta);
     (1+eta)*(1-zeta);
    -(1+eta)*(1-zeta);
    -(1-eta)*(1+zeta);
     (1-eta)*(1+zeta);
     (1+eta)*(1+zeta);
    -(1+eta)*(1+zeta);
];

dN_deta = 1/8 * [
    -(1-xi)*(1-zeta);
    -(1+xi)*(1-zeta);
     (1+xi)*(1-zeta);
     (1-xi)*(1-zeta);
    -(1-xi)*(1+zeta);
    -(1+xi)*(1+zeta);
     (1+xi)*(1+zeta);
     (1-xi)*(1+zeta);
];

dN_dzeta = 1/8 * [
    -(1-xi)*(1-eta);
    -(1+xi)*(1-eta);
    -(1+xi)*(1+eta);
    -(1-xi)*(1+eta);
     (1-xi)*(1-eta);
     (1+xi)*(1-eta);
     (1+xi)*(1+eta);
     (1-xi)*(1+eta);
];
end