function all_ip_map = project_tie_integration_points2(model, edata, step_idx)
% Project integration points from each slave (primary) surface face onto the master (secondary)
% mesh following FEBio's normal projection and fallback logic.
%
% Inputs:
%   model      - structure containing tie contact info
%   edata      - structure with node coordinates per step
%   step_idx   - simulation step to use
%
% Output:
%   all_ip_map - cell array {nface} where each cell is [nip x 1] struct array of projections

% -----------------
% DATA EXTRACTION
% -----------------
nodes = edata.steps{step_idx}.results.ecoords;            % [nnode x 3] (x, y, z)
primary_faces = model.tie_contact.primary_nodes;          % [prim x 6] (sec_face_id, n1-n4, hex_elem_id)
secondary_elems = model.tie_contact.secondary_elem;       % [nsec x 6] (sec_face_id, n1-n4, hex_elem_id)
nface = size(primary_faces, 1);

all_ip_map = cell(nface, 1);                              % output per face

face_centers = model.tie_contact.secondary_face_centers;  % [nsec x 4], already in the struct
Mdl = createns(face_centers(:,2:4), 'NSMethod', 'kdtree');% kd-tree to the center face
search_radius = 0.4;                                      % adjust with the typical mesh size (ex: 0.4x avg size)

for iface = 1:nface
    % --- Get coordinates of current primary face ---
    this_face_nodes = primary_faces(iface, 2:5);
    quad_xyz = nodes(this_face_nodes, :);         % [4 x 3]
    [ip_coords, ip_rs, ip_weights] = quad4_integration_points(quad_xyz);  % integration point info
    nip = length(ip_weights);


    ip_map(nip,1) = struct( ...
                        'prim_face',   [], ...
                        'prim_nodes',  [], ...
                        'prim_rs',     [], ...
                        'prim_xyz',    [], ...
                        'weight',      [], ...
                        'sec_face',    [], ...
                        'sec_elem_id', [], ...
                        'sec_rs',      [], ...
                        'proj_xyz',    [], ...
                        'dist',        [], ...
                        'gap',         [], ...
                        'method',      '', ...
                        'is_inside',   false, ...
                        'extrapolated', false);

    for ip = 1:nip
        % --- Integration point location and shape parameters ---
        P0 = ip_coords(ip, :);           % global xyz of IP
        rs_ip = ip_rs(ip, :);            % local (r,s) in primary face
        w_ip  = ip_weights(ip);          % integration weight

        % --- Compute primary surface normal at IP ---
        n_slave = quad4_normal_at_rs(quad_xyz, rs_ip);

       % -- kd-tree: search faces neart to the primary surface --
        [idx, dist] = rangesearch(Mdl, P0, search_radius);
        candidates = idx{1};  % index of the closest secondary faces

        % Initialize best match tracking
        min_dist_ray = Inf;
        min_t_ray = Inf;
        found_ray = false;

        % Track best (for fallback)
        min_dist_fallback = Inf;
        sec_face_fallback = [];
        rs_master_fallback = NaN(1,2);
        proj_xyz_fallback = NaN(1,3);

        % Variables for ray hit
        sec_face_ray = [];
        rs_master_ray = NaN(1,2);
        proj_xyz_ray = NaN(1,3);
        sec_elem_id_ray = NaN;

        % --- Search all secondary (master) faces ---
        for sec_i = candidates
            sec_nodes = secondary_elems(sec_i,2:5);   % node IDs (quad)
            sec_xyz = nodes(sec_nodes, :);            % coords [4 x 3]
            sec_elem_id = secondary_elems(sec_i,6);   % element id

            % ---- Ray-triangle intersection: split quad into two triangles ----
            % Triangle 1: nodes 1,2,4
            tri1 = sec_xyz([1,2,4], :);
            [hit1, t1, bary1] = ray_triangle_intersect(P0, n_slave, tri1);
            if hit1 && (t1 < min_t_ray) % Only hits in normal direction and closest
                rs_tri1 = bary1(2:3); % use second and third barycentric coords

                % Initial estimate of parametric coordinates in the quad,
                % based on triangle barycentric values (assume rs_tri from ray_triangle_intersect)
                rp = -1.0 + 2.0 * rs_tri1(1); % maps triangle (0,1) to quad (-1,+1)
                sp = -1.0 + 2.0 * rs_tri1(2);

                g = t1; % distance along the ray to intersection
                l1_guess = rp;
                l2_guess = sp;
                t_guess  = g;

                % Refine projection with Newton's method (more accurate intersection on quad)
                [rs_newton, t_newton, is_inside_newton] = newton_quad_projection(P0, n_slave, sec_xyz, l1_guess, l2_guess, t_guess);

                N_newton = quad_shape_functions(rs_newton(1), rs_newton(2)); 
                proj_xyz = N_newton * sec_xyz; % refined projected point on the quad surface

                rs = rs_newton;                % refined (r,s) coordinates
                is_inside = is_inside_newton;  % refined is-inside logical

                %[rs, is_inside] = get_quad_rs(proj_xyz, sec_xyz, sec_nodes); % Find (r,s) in quad
                if is_inside
                    found_ray = true;
                    min_t_ray = t_newton;
                    min_dist_ray = t_newton;
                    sec_face_ray = sec_nodes;
                    rs_master_ray = rs;
                    proj_xyz_ray = proj_xyz;
                    sec_elem_id_ray = sec_elem_id;
                end
            end
            % Triangle 2: nodes 3,4,2
            tri2 = sec_xyz([3,4,2], :);
            [hit2, t2, bary2] = ray_triangle_intersect(P0, n_slave, tri2);
            if hit2 && (t2 < min_t_ray)
                rs_tri2 = bary2(2:3);
                % Initial guesses for quad parametric coords from triangle 2 intersection
                rp = 1.0 - 2.0 * rs_tri2(1);  
                sp = 1.0 - 2.0 * rs_tri2(2);
                g  = t2;
            
                % Newton refinement
                [rs_newton, t_newton, is_inside_newton] = newton_quad_projection(P0, n_slave, sec_xyz, rp, sp, g);
            
                N_newton = quad_shape_functions(rs_newton(1), rs_newton(2));
                proj_xyz = N_newton * sec_xyz;
                % (Or: proj_xyz = P0 + t_newton * n_slave;)
            
                rs = rs_newton;
                is_inside = is_inside_newton;
                %[rs, is_inside] = get_quad_rs(proj_xyz, sec_xyz, sec_nodes);
                if is_inside
                    found_ray = true;
                    min_t_ray = t_newton;
                    min_dist_ray = t_newton;
                    sec_face_ray = sec_nodes;
                    rs_master_ray = rs;
                    proj_xyz_ray = proj_xyz;
                    sec_elem_id_ray = sec_elem_id;
                end
            end

            % --- Fallback: closest point on quad face ---
            [rs_cand, ~] = get_quad_rs(P0, sec_xyz, sec_nodes); % Newton's method
            N_cand = quad_shape_functions(rs_cand(1), rs_cand(2));
            proj_xyz_cand = N_cand * sec_xyz;
            dist_cand = norm(P0 - proj_xyz_cand);

            if dist_cand < min_dist_fallback
                min_dist_fallback = dist_cand;
                sec_face_fallback = sec_nodes;
                rs_master_fallback = rs_cand;
                proj_xyz_fallback = proj_xyz_cand;
            end
        end

        % === Final mapping selection (ray hit preferred over fallback) ===
        if found_ray
            method = 'ray';
            sec_face = sec_face_ray;
            rs_master = rs_master_ray;
            proj_xyz = proj_xyz_ray;
            gap = proj_xyz - P0;
            sec_elem_id = sec_elem_id_ray;
            is_inside = all(abs(rs_master)<=1+1e-8);
            extrapolated = false;
        else
            method = 'closest';
            sec_face = sec_face_fallback;
            rs_master = rs_master_fallback;
            proj_xyz = proj_xyz_fallback;
            gap = proj_xyz - P0;
            % Track if closest point is inside master parametric domain
            is_inside = all(abs(rs_master)<=1+1e-8);
            extrapolated = ~is_inside;
            sec_elem_id = NaN;
        end

        % Fill the map structure
        ip_map(ip).prim_face   = iface;           % index of primary face
        ip_map(ip).prim_nodes  = this_face_nodes; % node IDs (slave)
        ip_map(ip).prim_rs     = rs_ip;           % (r,s) on slave quad
        ip_map(ip).prim_xyz    = P0;              % xyz location of IP
        ip_map(ip).weight      = w_ip;            % integration weight

        ip_map(ip).sec_face    = sec_face;        % master quad node IDs
        ip_map(ip).sec_elem_id = sec_elem_id;     % master element id
        ip_map(ip).sec_rs      = rs_master;       % (r,s) on master
        ip_map(ip).proj_xyz    = proj_xyz;        % xyz projected onto master
        ip_map(ip).dist        = norm(gap);       % distance/slip magnitude
        ip_map(ip).gap         = gap;             % full vector gap

        ip_map(ip).method      = method;          % 'ray' (preferred) or 'closest' (fallback)
        ip_map(ip).is_inside   = is_inside;       % projected point is inside master parametric quad
        ip_map(ip).extrapolated= extrapolated;    % true if fallback is outside quad

    end
    all_ip_map{iface} = ip_map;
end
end




%%%%%%%%%%%%%%%%%%%%%%
% -- Helper Section --
%%%%%%%%%%%%%%%%%%%%%%

% Compute slave surface normal at (r,s)
function n = quad4_normal_at_rs(quad_xyz, rs)
    r = rs(1); s = rs(2);
    dN_dr = 0.25 * [-(1-s), (1-s), (1+s), -(1+s)];
    dN_ds = 0.25 * [-(1-r), -(1+r), (1+r), (1-r)];
    dx_dr = dN_dr * quad_xyz;
    dx_ds = dN_ds * quad_xyz;
    n = cross(dx_dr, dx_ds);
    n = n / norm(n);
end

function [intersect, t, bary] = ray_triangle_intersect(orig, dir, tri)
% Mimics FEBio C++ routine (metric tensor, dual basis, tolerance)
% Inputs:
%   orig: 1x3 ray origin
%   dir: 1x3 ray direction
%   tri: 3x3 triangle vertices, [v0; v1; v2]
% Outputs:
%   intersect: true/false, intersection found inside triangle (with eps tolerance)
%   t: distance along ray to intersection
%   bary: barycentric coordinates [b0, b1, b2]

    eps = 1e-4; % Tolerance, like FEBio's 'eps'

    % Edge vectors
    e0 = tri(2,:) - tri(1,:);
    e1 = tri(3,:) - tri(1,:);

    % Triangle normal (normalized)
    m = cross(e0, e1);
    m_norm = norm(m);
    if m_norm == 0
        intersect = false; t = NaN; bary = NaN(1,3);
        return;
    end
    m = m / m_norm;

    % Dot product of ray direction and triangle normal
    d = dot(dir, m);

    if d < 0 % only consider orientation "against" ray dir
        % Distance from ray origin to plane of triangle
        t = dot(m, tri(1,:) - orig) / d;

        % Intersection point
        q = orig + dir * t;

        % Metric tensor
        G = [dot(e0,e0), dot(e0,e1); dot(e1,e0), dot(e1,e1)];

        % Inverse metric tensor
        Gi = inv(G);

        % Dual basis
        E0 = Gi(1,1)*e0 + Gi(1,2)*e1;
        E1 = Gi(2,1)*e0 + Gi(2,2)*e1;

        % Parametric (local) coords
        rs0 = dot(E0, q-tri(1,:));
        rs1 = dot(E1, q-tri(1,:));

        % Barycentric coordinates as [v0, v1, v2]
        bary = [1 - rs0 - rs1, rs0, rs1];

        % Check if intersection point is within triangle (with tolerance)
        if (rs0 >= -eps) && (rs1 >= -eps) && (rs0 + rs1 <= 1+eps) 
            intersect = true;
            return;
        end
    end

    % If we reach here: no intersection
    intersect = false;
    t = NaN;
    bary = NaN(1,3);
end


% Find local (r,s) coordinates in QUAD4 for point Q (Newton iteration)
function [rs, is_inside] = get_quad_rs(Q, quad_xyz, ~)
    maxit = 350; tol = 1e-4;
    rs = [0,0];
    for it=1:maxit
        N = quad_shape_functions(rs(1),rs(2));
        Q_est = N*quad_xyz;
        F = Q_est - Q;
        if norm(F)<tol
            break;
        end
        dr = 1e-6;
        N1 = quad_shape_functions(rs(1)+dr, rs(2));
        Q1 = N1*quad_xyz;
        J1 = (Q1-Q_est)/dr;
        N2 = quad_shape_functions(rs(1), rs(2)+dr);
        Q2 = N2*quad_xyz;
        J2 = (Q2-Q_est)/dr;
        J = [J1(:), J2(:)];
        delta = J\F(:);
        rs = rs - delta';
    end
    is_inside = all(abs(rs)<=1+1e-4);
end

% QUAD4 shape functions
function N = quad_shape_functions(r,s)
    N = 0.25 * [(1-r)*(1-s), (1+r)*(1-s), (1+r)*(1+s), (1-r)*(1+s)];
end

% Generate QUAD4 integration points (2x2 Gauss rule)
function [ip_coords, ip_rs, ip_weights] = quad4_integration_points(quad_xyz)

    gauss_r = [-1, 1, 1,-1] * (1/sqrt(3));
    gauss_s = [-1,-1, 1, 1] * (1/sqrt(3));  % order: (ll, lr, ur, ul)

    gw = [1, 1, 1, 1]; % for 2x2 quadrature, all weights are 1
    
    ip_coords = zeros(4,3);
    ip_rs     = zeros(4,2); 
    ip_weights= ones(4,1); % all weights 1

    for ip = 1:4
        r = gauss_r(ip);
        s = gauss_s(ip);
        N = 0.25 * [(1 - r)*(1 - s), (1 + r)*(1 - s), (1 + r)*(1 + s), (1 - r)*(1 + s)];
        ip_coords(ip,:) = N * quad_xyz;
        ip_rs(ip,:)     = [r, s];
        ip_weights(ip)  = gw(ip);  % but all 1 for 2x2 scheme
    end
end


function [rs, t, is_inside] = newton_quad_projection(P0, dir, quad_nodes, l1_guess, l2_guess, t_guess)
% quad_nodes: 4x3 array, nodes in global 3D
% P0: ray origin, dir: direction
% l1_guess/l2_guess: initial guess for (r,s) [-1,+1]
% t_guess: initial guess for intersection time

    maxn = 5;
    eps_NR = 1e-7;
    l1 = l1_guess;
    l2 = l2_guess;
    l3 = t_guess;
    is_inside = false;

    for nn = 1:maxn
        % Shape functions and derivatives
        H  = 0.25 * [ (1-l1)*(1-l2); (1+l1)*(1-l2); (1+l1)*(1+l2); (1-l1)*(1+l2) ];
        dHdl1 = 0.25 * [ -(1-l2); (1-l2); (1+l2); -(1+l2)];
        dHdl2 = 0.25 * [ -(1-l1); -(1+l1); (1+l1); (1-l1)];
        
        % Quad position: sum_i H_i * node_i
        quad_pt = H' * quad_nodes;
        dquad_dl1 = dHdl1' * quad_nodes;
        dquad_dl2 = dHdl2' * quad_nodes;
        
        % Residual
        F = P0 + dir * l3 - quad_pt; % 1x3
        
        % Jacobian
        A = [ -dquad_dl1' -dquad_dl2' dir(:) ]; % 3x3
        
        % Solve update
        dx = -A \ F'; % column
        
        l1 = l1 + dx(1);
        l2 = l2 + dx(2);
        l3 = l3 + dx(3);
        
        if norm(dx) < eps_NR
            break;
        end
    end
    rs = [l1, l2];
    t = l3;
    % Inside quad
    if abs(l1) <= 1 + 1e-4 && abs(l2) <= 1 + 1e-4
        is_inside = true;
    end
end