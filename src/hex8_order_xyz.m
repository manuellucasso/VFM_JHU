function pass = hex8_order_xyz_signs(node_xyz,elem_id)
xi_eta_zeta = [
    -1 -1 -1;
     1 -1 -1;
     1  1 -1;
    -1  1 -1;
    -1 -1  1;
     1 -1  1;
     1  1  1;
    -1  1  1;
];

hex8_shape_functions = @(xi_eta_zeta) (1/8) * [ ...
    (1-xi_eta_zeta(1))*(1-xi_eta_zeta(2))*(1-xi_eta_zeta(3));
    (1+xi_eta_zeta(1))*(1-xi_eta_zeta(2))*(1-xi_eta_zeta(3));
    (1+xi_eta_zeta(1))*(1+xi_eta_zeta(2))*(1-xi_eta_zeta(3));
    (1-xi_eta_zeta(1))*(1+xi_eta_zeta(2))*(1-xi_eta_zeta(3));
    (1-xi_eta_zeta(1))*(1-xi_eta_zeta(2))*(1+xi_eta_zeta(3));
    (1+xi_eta_zeta(1))*(1-xi_eta_zeta(2))*(1+xi_eta_zeta(3));
    (1+xi_eta_zeta(1))*(1+xi_eta_zeta(2))*(1+xi_eta_zeta(3));
    (1-xi_eta_zeta(1))*(1+xi_eta_zeta(2))*(1+xi_eta_zeta(3))
]';

tolerance = 1e-10;
pass = true;

for k = 1:8
    N = hex8_shape_functions(xi_eta_zeta(k,:));     % Shape function weights at parametric node k
    node_interp = N * node_xyz;                    % Interpolate position at node k
    node_actual = node_xyz(k,:);
    err = norm(node_interp - node_actual);

    if err > tolerance
        pass = false;
    end
end

if pass
else
    disp('Node order DOES NOT match canonical HEX8 parametric mapping in element %d',elem_id);
end
end