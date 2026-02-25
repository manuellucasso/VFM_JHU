function [ops,A,B] = create_ops(target_rows, target_cols, source_rows, source_cols, weights,mat_size)
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

    [A, B] = build_AB(mat_size, ops);
end

function [A, B] = build_AB(mat_size, ops)
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

    n_rows = mat_size(1);
    n_cols = mat_size(2);
    n_ops  = numel(ops);

    % Start with identity matrix for B (keeps original values)
    B = eye(n_rows, n_rows);
    B = kron(B, ones(1, n_cols)); % expand to each column if needed
    B = reshape(B, n_rows, n_cols); % make sure B is the correct size

    % Start with zeros for A (will only set source–>target values)
    A = zeros(n_rows, n_cols);

    % Loop through all operations
    for k = 1:n_ops
        tgt = ops(k).target;     % target row index
        src = ops(k).source;     % source row index
        ll  = ops(k).cols;       % vector of columns to set

        % For each specified column:
        for c = ll
            % Set corresponding element in A: take from the source row
            A(tgt, c) = src;

            % Set B so the original value at that [tgt,c] is ignored
            B(tgt, c) = 0;
        end
    end

end