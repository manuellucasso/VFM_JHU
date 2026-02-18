function [model, edata] = preliminary_reading(mydir, mymodel, myexpdata,edata,model)
% Reads a FEBio model and corresponding experimental (or simulation) data.
% Returns:
%   model - finite element model struct (geometry, connectivity, etc.)
%   edata - struct containing time steps, extracted results per step



% --- Load FEBio model (should contain elements, nodes, surfaces info) ---
if (~exist('model', 'var')) || isempty(fieldnames(model))
    % 'model' does not exist OR is an empty struct (has no fields)
    model = read_febio_model_test(mydir, mymodel);
else
    disp('model already exists and is not empty. Not reading.');
end

nel = length(model.elements);        % Number of elements
nnd = length(model.nodes);           % Number of nodes
nels = length(model.surfacesp);      % Number of surfaces (possibly unused here)

% --- Open the experimental data file ---
fullpath_expdata = fullfile(mydir, myexpdata);
fid = fopen(fullpath_expdata, 'r');
line = fgetl(fid);                   % Read first line

% --- Initialize variables to store step info and times ---
steps = {};        % Will become cell array of structs for each time step
times = [];        % Times (float) at each step
step_idx = 0;      % Current step index
last_time_value = NaN; % Tracks last time value to avoid duplicates

% --- Read file line by line ---
while ischar(line)
    % --- If the line marks a new time step, log the time and start a new step ---
    if strncmp(line, 'Time =', 6)
        this_time = sscanf(line, 'Time = %f');
        % Only increment if we see a new time value
        if isempty(times) || this_time ~= last_time_value
            step_idx = step_idx + 1;
            steps{step_idx} = struct();     % New step struct
            times(step_idx) = this_time;    % Store time value
            last_time_value = this_time;    % Update last time
        end
        line = fgetl(fid);  % Move to next line
        continue
    end

    % --- Node data block: coordinates and displacements ---
    if strcmp(line, 'Data = x;y;z;ux;uy;uz')
        % Prepare to read per-node data: node ID, coords, displacements
        nodedat = zeros(nnd,7); % Columns: id, x, y, z, ux, uy, uz
        for i = 1:nnd
            nodedat(i,:) = fscanf(fid, '%i %f %f %f %f %f %f\n', [1, 7]);
        end
        % Store reference coordinates and displacements in this step
        steps{step_idx}.results.ecoords = nodedat(:,2:4);  % (x,y,z)
        steps{step_idx}.results.edisp   = nodedat(:,5:7);  % (ux,uy,uz)

    end

    % --- Move to next line for next iteration ---
    line = fgetl(fid);
end

% --- Close file after reading ---
fclose(fid);

% --- Package results into output struct (edata) ---
%   edata.steps: cell array, each cell is a struct for a time step.
%   edata.times: vector of time values for each step.
edata.steps = steps;
edata.times = times;


num_steps = [1];


% Só inicializa se ainda não existe
if ~isfield(edata, 'all_ip_map')
    edata.all_ip_map = cell(2, 1);
end


% In case there is no contact, skip the project
for i = 1:length(num_steps)
    if isfield(model, 'tie_contact') && ~isempty(model.tie_contact)
        if isempty(edata.all_ip_map{i})  % Só roda se está vazio!
            edata.all_ip_map{i} = project_tie_integration_points2(model, edata, num_steps(i));
        else
            disp(['Mapping for step ', num2str(num_steps(i)), ' already exists. Skipping.']);
        end
    else
        disp('Skipping mapping: tie_contact.contact is missing or empty.');
    end
end

% Retrieve the all_ip_map for that step
all_ip_map = edata.all_ip_map{1};


end