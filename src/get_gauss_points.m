function [gauss_points, weights] = get_gauss_points(dimension, order)
% GET_GAUSS_POINTS Unified function for Gauss points in 1D, 2D, and 3D
%   dimension: 1, 2, or 3 (for line, quadrilateral, or hexahedral elements)
%   order: integration order (1, 2, 3, etc.)

    % 1D Gauss points and weights
    switch order
        case 1
            xi_1d = 0;
            w_1d = 2;
        case 2
            xi_1d = [-1/sqrt(3), 1/sqrt(3)];
            w_1d = [1, 1];
        case 3
            xi_1d = [-sqrt(3/5), 0, sqrt(3/5)];
            w_1d = [5/9, 8/9, 5/9];
        otherwise
            error('Unsupported Gauss order');
    end
    
    % Generate points and weights based on dimension
    switch dimension
        case 1  % 1D - Line elements
            gauss_points = xi_1d(:);
            weights = w_1d(:);
            
        case 2  % 2D - Quadrilateral elements
            [xi_grid, eta_grid] = meshgrid(xi_1d, xi_1d);
            gauss_points = [xi_grid(:), eta_grid(:)];
            
            [wxi, weta] = meshgrid(w_1d, w_1d);
            weights = wxi(:) .* weta(:);
            
        case 3  % 3D - Hexahedral elements
            [xi_grid, eta_grid, zeta_grid] = meshgrid(xi_1d, xi_1d, xi_1d);
            gauss_points = [xi_grid(:), eta_grid(:), zeta_grid(:)];
            
            [wxi, weta, wzeta] = meshgrid(w_1d, w_1d, w_1d);
            weights = wxi(:) .* weta(:) .* wzeta(:);
            
        otherwise
            error('Dimension must be 1, 2, or 3');
    end
end