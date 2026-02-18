function [datau,dataF]=simulate_febio_uniform(mydir,mymodel,matparam,nnd,nel,param_ind,changing_matrix)
global last_time
 
%dimension storage for simulation results
%datau=storage for displacement (u)
datau = zeros(nnd,7);

%dataF=storage for deformation gradient (F) and discplacement gradient (G)
dataF = zeros(nel,19);

%simulate data with target properties. It reads in template.
%nummat=Number of materials to modify
nummat = size(changing_matrix,2);

% Read file as before
fid = fopen([mydir '/' mymodel],'r');
data = textscan(fid, '%s', 'Delimiter', '\n', 'CollectOutput', true);
fclose(fid);
data = data{1}; % flatten cell array
Nline = numel(data);



for matn = 1:nummat
    mat_id = changing_matrix(1,param_ind);
    material_start_line = 0;
    material_end_line = 0;
    % Removed is_unique and count_here, no longer needed

    % --- Identify block ---
    for i = 1:Nline
        if contains(data{i}, sprintf('material id="%d"', mat_id))
            material_start_line = i;
        end
        if material_start_line > 0 && contains(data{i}, '</material>')
            material_end_line = i;
            break;
        end
    end

% --- If block found ---
    if material_start_line > 0 && material_end_line > 0
        % --- Find material type for this block ---
        blocktype = '';
        for searchline = material_start_line:material_end_line
            if contains(data{searchline}, 'neo-Hookean')
                blocktype = 'neo-Hookean';
                break;
            elseif contains(data{searchline}, 'type="coupled Mooney-Rivlin')
                blocktype = 'Mooney-Rivlin';
                break;
            elseif contains(data{searchline}, 'type="coupled trans-iso Mooney-Rivlin">')
                blocktype = 'coupled trans-iso Mooney-Rivlin';
                break;
            elseif contains(data{searchline}, 'type="HGO unconstrained"')
                blocktype = 'HGO unconstrained';
                break;
            end
            % Add more elseif for other types if needed
        end

        for j = material_start_line:material_end_line
            switch blocktype
                case 'neo-Hookean'
                    % --- Neo-Hookean: change both <E> and <v> ---
                    if contains(data{j}, '<E>')
                        data{j} = sprintf('<E>%.8g</E>', matparam(param_ind, 1)); % E
                    end
                    if contains(data{j}, '<v>')
                        data{j} = sprintf('<v>%.8g</v>', matparam(param_ind, 3)); % v
                    end

                case 'Mooney-Rivlin'
                    % --- Mooney-Rivlin: change <c1>, <k> ---
                    if contains(data{j}, '<c1>')
                        data{j} = sprintf('<c1>%.8g</c1>', matparam(param_ind, 1));
                    end
                    if contains(data{j}, '<k>')
                        data{j} = sprintf('<k>%.8g</k>', matparam(param_ind, 3));
                    end

                case 'coupled trans-iso Mooney-Rivlin'
                    % --- Change all parameters for trans-iso Mooney-Rivlin ---
                    if contains(data{j}, '<c1>')
                        data{j} = sprintf('<c1>%.8g</c1>', matparam(param_ind, 1));
                    elseif contains(data{j}, '<c2>')
                        data{j} = sprintf('<c2>%.8g</c2>', matparam(param_ind, 2));
                    elseif contains(data{j}, '<k>')
                        data{j} = sprintf('<k>%.8g</k>', matparam(param_ind, 3));
                    elseif contains(data{j}, '<c3>')
                        data{j} = sprintf('<c3>%.8g</c3>', matparam(param_ind, 4));
                    elseif contains(data{j}, '<c4>')
                        data{j} = sprintf('<c4>%.8g</c4>', matparam(param_ind, 5));
                    elseif contains(data{j}, '<c5>')
                        data{j} = sprintf('<c5>%.8g</c5>', matparam(param_ind, 6));
                    elseif contains(data{j}, '<lam_max>')
                        data{j} = sprintf('<lam_max>%.8g</lam_max>', matparam(param_ind, 7));
                    end

                case 'HGO unconstrained'
                    % --- HGO unconstrained: update <c>, <k1>, <k>, <k2>, <kappa>, <gamma> ---
                    if contains(data{j}, '<c>')
                        data{j} = sprintf('<c>%.8g</c>', matparam(param_ind, 1));       % c
                    elseif contains(data{j}, '<k1>')
                        data{j} = sprintf('<k1>%.8g</k1>', matparam(param_ind, 2));     % k1
                    elseif contains(data{j}, '<k>')
                        data{j} = sprintf('<k>%.8g</k>', matparam(param_ind, 3));       % k
                    elseif contains(data{j}, '<k2>')
                        data{j} = sprintf('<k2>%.8g</k2>', matparam(param_ind, 4));     % k2
                    elseif contains(data{j}, '<kappa>')
                        data{j} = sprintf('<kappa>%.8g</kappa>', matparam(param_ind, 5)); % kappa
                    elseif contains(data{j}, '<gamma>')
                        data{j} = sprintf('<gamma>%.8g</gamma>', matparam(param_ind, 6)); % gamma
                    end
                % Add more cases here for other material types as needed
            end
        end
    end
end

fid = fopen([mydir '/' 'modeltilde.feb'], 'w');
for I = 1:Nline
    fprintf(fid, '%s\n', char(data{I}));
end
fclose(fid);

%run febio simulations with new material properties.
%status 0=Success (FEBio ran without errors)
%status 1=Generic error (e.g., Convergence problem)
%status -1=Command could not be executed (e.g., path incorrect)
febio_bin = '/home/msampai4/FEBio/febio_src/build/bin/febio4';
unique_model_name = 'modeltilde.feb';
unique_log_name = 'modeltilde.log'; 
unique_output_name = 'output.txt'; % <-- unique filename
febio_cmd = sprintf('"%s" -i %s > %s 2>&1', ...
    febio_bin, unique_model_name, unique_output_name);

[status, cmdout] = system(febio_cmd);

%% Check for convergence
% Check for failure conditions
if status ~= 0 || contains(cmdout, 'E R R O R   T E R M I N A T I O N')
    error_flag = true;
    disp('FEBio terminated with error after completing timing info');
    
    % Match original nodedat dimensions but fill with NaN
    datau = nan(nnd, 7);    
    
    % Match original elemdat dimensions but fill with NaN
    dataF = nan(nel, 19);   
    
    %leave the function after the error
    return
else
    error_flag = false;
end


%% Parse FEBio's log file o extract nodal displacement and deformation gradients

% Open file for reading
fres = fopen('modeltilde.log','r');
if fres == -1
    error('Could not open log file.');
end

% Construct string to search for Time
time_string = sprintf('Time = %.6g', last_time);

% -- Search for the Time header --
s = fgetl(fres);
found_time = false;
while ischar(s)
    % Robust search: also matches slightly formatted times (e.g., 1.5 == 1.500)
    if contains(s, 'Time =')
        tmp_time = sscanf(s, 'Time = %f');
        if ~isempty(tmp_time) && abs(tmp_time - last_time) < 1e-8
            found_time = true;
            break;
        end
    end
    s = fgetl(fres);
end

if ~found_time
    fclose(fres);
    error('Requested time %g not found in log file.', last_time);
end

% -- Search for Displacement Data header --
displ_header = 'Data = x;y;z;ux;uy;uz';
while ischar(s)
    if strcmp(s, displ_header)
        break;
    end
    s = fgetl(fres);
end
if ~ischar(s)
    fclose(fres);
    error('Displacement data header not found.');
end

% -- Parse displacement data --
datau = zeros(nnd, 7);
for idk = 1:nnd
    vals = fscanf(fres,'%i %f %f %f %f %f %f\n',[1,7]);
    if isempty(vals) || numel(vals) ~= 7
        fclose(fres);
        error('Displacement data incomplete at node %d.', idk);
    end
    datau(idk, :) = vals;
end


fclose(fres);

