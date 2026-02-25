function ops = create_ops(target_rows, target_cols, source_rows, source_cols, weights)
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
end