function fiber_table = read_febio_fibers_elementset(xDoc, elem_set_name)
% Reads fiber directions for specified element set from XML document

elBlocks = xDoc.getElementsByTagName('Elements');
part_elem_ids = [];
for b = 0:elBlocks.getLength-1
    elemBlock = elBlocks.item(b);
    if strcmp(char(elemBlock.getAttribute('name')), elem_set_name)
        elems = elemBlock.getElementsByTagName('elem');
        Nel = elems.getLength;
        part_elem_ids = zeros(Nel, 1);
        for e = 0:Nel-1
            elem = elems.item(e);
            eid = str2double(elem.getAttribute('id'));
            part_elem_ids(e+1) = eid;
        end
        break
    end
end

if isempty(part_elem_ids)
    error('Element set "%s" not found in <Elements>', elem_set_name);
end

% Find corresponding ElementData node
eldata_nodes = xDoc.getElementsByTagName('ElementData');
fiber_eldata = [];
for k = 0:eldata_nodes.getLength-1
    nd = eldata_nodes.item(k);
    if strcmp(char(nd.getAttribute('name')), 'fiber') && ...
       strcmp(char(nd.getAttribute('elem_set')), elem_set_name)
        fiber_eldata = nd;
        break
    end
end

if isempty(fiber_eldata)
    error('Fiber ElementData for elem_set="%s" not found', elem_set_name);
end

fiber_elems = fiber_eldata.getElementsByTagName('e');
Nfib = fiber_elems.getLength;
fiber_table = zeros(Nfib, 4); % [element_id, fiber_x, fiber_y, fiber_z]
for i = 0:Nfib-1
    e_node = fiber_elems.item(i);
    lid = str2double(e_node.getAttribute('lid')); % starts at 1
    vec_str = char(e_node.getTextContent);
    fiber_vec = sscanf(vec_str, '%f,%f,%f');
    if numel(fiber_vec) ~= 3
        fiber_vec = sscanf(vec_str, '%f %f %f'); % fallback
    end
    elem_id = part_elem_ids(lid);
    fiber_table(i+1,:) = [elem_id, fiber_vec(:)'];
end

end