%--------------------------------------------------------------------------
% Main Script: Virtual Field Method Parameter Optimization and Cost Function Analysis
%--------------------------------------------------------------------------
% This script conducts a parameter sweep for Neo-Hookean material model fitting using
% Internal Virtual Work and External Virtual Work techniques and stores results.
% It uses Latin Hypercube Sampling, calculates cost function, performs post-processing, 
% and writes results to CSV.
%--------------------------------------------------------------------------

clear all;
clc;
fclose all;

% Declare global variables for tracking across optimization runs
global totalRunCount IVW_vector EVW_vector last_time TieVW_vector

%% --------------------------------------------------------------------------
% Directories and Data Files
%--------------------------------------------------------------------------
% Folder containing your FE model and node/element data logs
mydir = 'C:\Users\manue\OneDrive - University of Ottawa\JHU\Research\Virtual Field Methods\File\VFM\Matlab_2Variables\Tensile Uniform Load\Two Cost Function\Presstressed_Model\Fpre_calc\';

% FEBio model file with mesh, boundary, initial/material configuration
mymodel = 'Not_matching.feb'; 

% Log file containing "experimental node data" and calculated element strains
myexpdata = 'Not_matching.log'; 

matFile = 'experiment_results_parsed.mat'; 

preFile = 'experimental_Fpre.mat';

%% --------------------------------------------------------------------------
% Material Parameter Bounds and Reference Pressure
%--------------------------------------------------------------------------

% Reference and applied pressure (p_ref is unused here) 
p_app = -0.003;          % Applied pressure
%p_app=[0.15;0.02];


% Last physical time at which prestress is calculated/frozen
prestress_time = 1.0;

% Last Time of simulation
last_time = 2;
%last_time = 1.5;

% Penalty Factor for the contact
eps = 1000;

%% --------------------------------------------------------------------------
% Gauss Order for Numerical Integration
%--------------------------------------------------------------------------
gauss_order = 2;       % Selects 2x2x2 Gauss quadrature for elements

%% --------------------------------------------------------------------------
% Latin Hypercube Sampling for Design of Experiments
%--------------------------------------------------------------------------
nSamples = 50;        % Number of random parameter samples

nMaterial = 8;       % Number of  different material models

% Parameter bounds [lower, upper] for c1 and c2 (shear moduli)
lb = [0.1, 0.2,0.1,0.1,0.2,0.2,0.2,0.2];     % Lower bounds for [c1, c2]
ub = [0.3, 0.5,0.3,0.3,0.5,0.5,0.5,0.5];     % Upper bounds for [c1, c2]
%lb = [0.05, 0.3];   % lower bounds
%ub = [0.3,  1.0];   % upper bounds


% Ground-truth parameters
gt = [0.25,0.36,0.3,0.25,0.3,0.3,0.3,0.3];
%gt = [0.1,0.5];

% Virtual Work Parameter Variation
virtual_work_parameters = [0.1729824472142898,0.254981153753818,0.204602575103828,...
    0.1511037795764970,0.256911776087602,0.393591407063561,0.445408049930159,0.288630749044160];

%virtual_work_parameters = [0.212829455,	0.554808007] ;   

% Generate Latin Hypercube samples (efficient space-filling design)
lhs_matrix = lhsdesign(nSamples, nMaterial);

% Scaling to c_i range
for i=1:nMaterial
    lhs_matrix(:,i) = lb(i) + (ub(i) - lb(i)) * lhs_matrix(:,i);
end

% Including the ground-truth for verification
lhs_matrix = [lhs_matrix;gt];

parameters = sortrows(lhs_matrix, 1);

parameters = [virtual_work_parameters;parameters];



%% --------------------------------------------------------------------------
% Allocate Storage for Results
%--------------------------------------------------------------------------
nruns_c1 = length(parameters(:,1));
nruns = nruns_c1;                                         % Total number of parameter combinations

IVW = zeros(1, nruns);                                    % Internal Virtual Work (per run)
EVW = zeros(1, nruns);                                    % External Virtual Work (per run)
TieVW = zeros(1, nruns);                                  % Tie Contact Virtual Work (per run)
cost_func = zeros(1, nruns);                              % Cost function: |IVW - EVW|^2 (per run)
model = struct();                                         % model data (nodes, elements, surfaces, etc.)
edata = struct();                                         % simulation data (experimental data or results)



totalRunCount = 2                                        % Run counter

%% --------------------------------------------------------------------------
% Read Model and Node/Element Data, Compute Cumulative Fpre
%  --------------------------------------------------------------------------
% Loads mesh and log file, including node coordinates, element connectivity, results at timesteps
if isfile(matFile)
    % If the results file exists, load the data
    loaded = load(matFile);           % Loads to a struct 'loaded'
    model = loaded.model;
    edata = loaded.edata;
    disp('Loaded model and edata from saved .mat file.');
else
    [model, edata] = preliminary_reading(mydir, mymodel, myexpdata,edata,model);
    save(matFile, 'model', 'edata');
    disp('Read experiment dataparsed and saved results.');
end



%% --------------------------------------------------------------------------
% Main Parameter Sweep Loop
%--------------------------------------------------------------------------
for i = 1:nruns
    for k = 1:nMaterial  

        % Assign material properties for plate and inclusion for this run
        matparam(k,1) = parameters(i,k);
        matparam(k,2) = parameters(i,k) * 10000;

    end

    % Compute cumulative Fpre up to prestress_time, and freeze beyond
    edata = accumulate_Fpre_from_edata(mydir, mymodel, gauss_order, prestress_time,matparam,edata);

    % Calculate virtual work quantities and cost function for current parameter set
    % IVW: Internal Virtual Work, EVW: External Virtual Work
    % cost_func: squared error of IVW vs. EVW (to be minimized)
    [IVW(totalRunCount), EVW(totalRunCount),TieVW(totalRunCount) ,cost_func(totalRunCount)] = ...
        calc_virtual_work_variation_integration2(mydir, mymodel, model, edata, matparam, p_app, gauss_order,eps);

    
    % Store material parameters for this run
    param_row = matparam(:,1)'
    if ~exist('paramRuns', 'var') || isempty(paramRuns)
        paramRuns = param_row;
    else
        paramRuns = [paramRuns; param_row];
    end
    % Increment run counter
    totalRunCount = totalRunCount + 1
    
end

%% --------------------------------------------------------------------------
%  Postprocessing: Find Minimum, Prepare Output, File Write, Visualization
%  --------------------------------------------------------------------------

% ----- Combine vectors into table for output -----

% Generate parameter column names automatically: c1, c2, ..., cn
paramHeaders = cell(1, nMaterial);
for k = 1:nMaterial
    paramHeaders{k} = ['c' num2str(k)];
end


% Other result columns to append
otherHeaders = {'EVW', 'IVW', 'TieVW', 'cost_func'};
otherVars = {transpose(EVW(2:end)), transpose(IVW(2:end)), transpose(TieVW(2:end)), transpose(cost_func(2:end))};

% Combine parameter columns and other columns into cell array
data_cells = cell(1, nMaterial + numel(otherVars));
for k = 1:nMaterial
    data_cells{k} = paramRuns(:, k);
end
for k = 1:numel(otherVars)
    data_cells{nMaterial + k} = otherVars{k};
end

% Combine all headers
allHeaders = [paramHeaders, otherHeaders];

% Create the final table with dynamic headers
data_table = table(data_cells{:}, 'VariableNames', allHeaders);

% Optionally, write to CSV
writetable(data_table, 'virtual_work_results.csv');

% ----- Plot: Cost Function vs Modulus c2 ------
c1= paramRuns(:,1);
c2= paramRuns(:,2);

figure;
loglog(c2(1:end), cost_func(2:end), 'b-o', 'LineWidth', 1.5, 'MarkerSize', 6);
grid on;
xlabel('Modulus c2', 'FontSize', 12);
ylabel('Cost Function abs(IVW-EVW)^2', 'FontSize', 12);
title('Cost Function vs Modulus (Log Scale)', 'FontSize', 14);

% ----- Plot: Cost Function vs Modulus c1 ------
figure;
loglog(c1(1:end), cost_func(2:end), 'b-o', 'LineWidth', 1.5, 'MarkerSize', 6);
grid on;
xlabel('Modulus c1', 'FontSize', 12);
ylabel('Cost Function abs(IVW-EVW)^2', 'FontSize', 12);
title('Cost Function vs Modulus (Log Scale)', 'FontSize', 14);

% ----- Prepare and filter data for 3D visualization (log scale) -----
valid = (c1 > 0) & (c2 > 0) & (cost_func > 0);
c1_valid = c1(1:end);
c2_valid = c2(1:end);
cost_func_valid = cost_func(2:end);

% ----- Delaunay triangulation for 3D surface plot -----
tri = delaunay(c1_valid, c2_valid);

figure;
trisurf(tri, c1_valid, c2_valid, cost_func_valid, ...
    'FaceColor', [0.2 0.6 0.8], ...
    'EdgeColor', 'none');
xlabel('c1'); ylabel('c2'); zlabel('cost\_func');
title('Cost Function vs. Shear Moduli (Log Scale on All Axes)');
set(gca, 'XScale', 'log', 'YScale', 'log', 'ZScale', 'log');    % Log scale axes
light('Position', [1 1 1]); lighting phong; material dull;
view(3); rotate3d on;