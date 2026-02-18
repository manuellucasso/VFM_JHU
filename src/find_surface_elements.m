function surfacesp_out = find_surface_elements(surfacesp, elements)
% Returns a vector with the same number of rows as surfacesp.
% Each entry gives the element ID that contains the nodes of that surface.

n_surfaces = size(surfacesp, 1);
surfacesp_out = zeros(n_surfaces, size(surfacesp, 2) + 1);  % initialize output

% Create a map from element node sets for faster comparison
n_elements = size(elements, 1);
element_node_sets = cell(n_elements, 1);
element_ids = elements(:, 1);  % first column is element ID

for i = 1:n_elements
    element_node_sets{i} = sort(elements(i, 2:end));  % sort node list for robust comparison
end

for s = 1:n_surfaces
    surf_id = surfacesp(s, 1);
    surf_nodes = sort(surfacesp(s, 2:end));  % remove surface ID, keep only nodes

    found = false;
    for e = 1:n_elements
        % if all surface nodes are in the element node list
        if all(ismember(surf_nodes, element_node_sets{e}))
            surfacesp_out(s, :) = [surfacesp(s, :) element_ids(e)];  % store element ID
            found = true;
            break;  % assume each surface belongs to one element
        end
    end

    if ~found
        warning('Surface ID %d does not belong to any element.', surf_id);
        surfacesp_out(s, :) = [surfacesp(s, :) NaN];
    end
end
end