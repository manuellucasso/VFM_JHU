function ops_matrix_struct = create_ops(target_rows, target_cols, source_rows, source_cols, weights,mat_size)
%CREATE_OPS_STRUCT Construct an array of operation structs for parameter mapping.
%
%   ops = CREATE_OPS_STRUCT(target_rows, target_cols, source_rows, source_cols, weights)
%
%   This function builds an array of structs where each struct describes
%   how a target parameter will be set as a linear combination of source parameters.
%
%   Inputs:
%       target_rows  - vector of target rows (one per operation)
%       target_cols  - vector of target columns (one per operation)
%       source_rows  - cell array, each cell contains vector of source row indices
%       source_cols  - cell array, each cell contains vector of source col indices
%       weights      - cell array, each cell contains vector of weights (same size as sources)
%
%   Output:
%       ops - array of struct with fields:
%             .target_row, .target_col, .source_rows, .source_cols, .weights

    % Number of operations to create
    n_ops = length(target_rows);

    % Preallocate the ops struct array
    ops = struct('target_row', {}, 'target_col', {}, ...
                 'source_rows', {}, 'source_cols', {}, 'weights', {});

    for k = 1:n_ops
        % Fill struct for each operation
        ops(k).target_row = target_rows(k);
        ops(k).target_col = target_cols(k);
        ops(k).source_rows = source_rows{k};
        ops(k).source_cols = source_cols{k};
        ops(k).weights = weights{k};
    end

    ops_matrix_struct = build_AB(mat_size, ops);
end

function ops_matrix_struct = build_AB(mat_size, ops)
%BUILD_AB Constructs selection matrices for parameter transfer.
%   [A,B] = BUILD_AB(mat_size, ops)
%
%   INPUTS:
%       mat_size - Size of the parameter matrix [rows, cols].
%       ops      - Array of structs, each with fields:
%                  .target   - Target row to receive values
%                  .source   - Source row to copy values from
%                  .cols     - Vector of columns (parameters) to modify
%
%   OUTPUTS:
%       A - Matrix (rows x cols), for the new values (from source rows)
%       B - Matrix (rows x cols), for original values (not modified by operations)

    % n_rows = number of materials
    n_rows = mat_size(1);
    
    % n_cols = number of paremeters per material
    n_cols = mat_size(2);
    
    % number of operation change
    n_ops  = numel(ops);

    % Start with identity matrix for B (keeps original values)
    B = eye(n_cols);

    % Start with zeros for A (will only set source–>target values)
    A = zeros(n_cols, n_cols);

    % Creating ops_matrix_struct with A and B matrix
    ops_matrix_struct = struct('A', {},'B', {},'tgt_row',{},'src_row',{});

    for i = 1:n_rows

        % Store A and B in the struct for this material (row)
        ops_matrix_struct(i).B = B;
        ops_matrix_struct(i).tgt_row = i;

    end

    % Loop through all operations
    for k = 1:n_ops
        tgt = ops(k).target_row;     % target row index
        src = ops(k).source_rows;       % source row index
        tgt_col  = ops(k).target_col;      % vector of columns to set
        src_col  = ops(k).source_cols;        % vector of source columns 
        weights  = ops(k).weights;               % weights of A matrix

        % Set corresponding source row
        ops_matrix_struct(tgt{1}).src_row = src;
        
        for i = 1:length(src)
        
            % Creating a matrix A for each source row    
            ops_matrix_struct(tgt{1}).A{i} = A;

            % For each specified column:
            for c = 1:length(tgt_col{1})          
                
                % Set corresponding element in A: take from the source row
                tgt_column_spc = tgt_col{1}(c);
                src_column_spc  = src_col(c);
    
                ops_matrix_struct(tgt{1}).A{i}(tgt_column_spc, src_column_spc) = 1*weights(src_column_spc);
                
                % Set B so the original value at that [tgt,c] is ignored
                ops_matrix_struct(tgt{1}).B(tgt_column_spc, src_column_spc) = 0;
            end
        end
    end

end