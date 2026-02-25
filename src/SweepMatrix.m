function [ground_truth_mat,matparam_sweep,matparam_complete] = SweepMatrix(model,changing_matrix,x,Normalizer,ops_matrix_struct)
    % Compute Cauchy stress for material in 'model' at index 'mat_idx'
    % F: 3x3 deformation gradient
    % model: struct from your parser (see above)
    % mat_idx: index of material to use
           
    ground_truth_mat = model.matprop;
    matparam_complete = ground_truth_mat;
    corresponding = cell2mat(changing_matrix(1,:));
    parameters = x .* Normalizer;
    matparam_sweep = zeros([length(corresponding),size(ground_truth_mat,2)]);


    %Creating material parameters vector to pass on to the next functions
    % Assign material properties for this run
    for k = 1:max(corresponding)  

        if ismember(k,corresponding)
            idx_int = find(corresponding == k, 1, 'first');
            idx_last = find(corresponding == k, 1, 'last');
            nParam = idx_last - idx_int + 1;
            
            if idx_int-idx_last == 0
                % The retina is not incompressible
                if k == 3
                    matparam(k,1) = parameters(idx_int);
                    matparam(k,2) = parameters(idx_int) * 20;
                else
                
                    matparam(k,1) = parameters(idx_int);
                    matparam(k,2) = parameters(idx_int) * 1000;
                end

            else
                for iParam = 1:nParam
                    matparam(k,iParam) = parameters(idx_int+iParam-1);
                end

            end
        end

    end

    for col = 1:size(changing_matrix, 2)
        % Determine material model string
        mat_idx = changing_matrix{1,col};
        mat_type = lower(strtrim(model.matmodel(mat_idx)));
        mat_change = changing_matrix{4,col};
        
        % Getting the order of the property for a given material type
        prop = mat_change2prop(mat_type, mat_change);

        % Changing the value in the material parameter sweep matrix
        col_param_prop = prop;
        row_param = changing_matrix{1,col}; 
        col_param = changing_matrix{2,col};  
        matparam_sweep(col,:) = ground_truth_mat(row_param,:);
        matparam_sweep(col,col_param_prop) = matparam(row_param,col_param);


        % Changing the value in the material parameter matparam_complete
        matparam_complete(mat_idx,col_param_prop) = matparam(row_param,col_param);
   
    end

   
    % Performing operations between parameters
    totalNparam =length(matparam_complete(1,:));
    %Changing the K column from the matparam_complete
    totalNmat = length(matparam_complete(:,1));

    for idx_Mat = 1:totalNmat
        tgt_row = ops_matrix_struct(idx_Mat).tgt_row;
        src_row = ops_matrix_struct(idx_Mat).src_row;
        A_matrix = ops_matrix_struct(idx_Mat).A;
        B_matrix = ops_matrix_struct(idx_Mat).B;

        src_mult = zeros(1,totalNparam);

        if isempty(A_matrix)
            continue
        end

        for n_src = 1:length(src_row)
        
           src_mult = src_mult + matparam_complete(src_row(n_src),:)* A_matrix{n_src}; 
        
        end

        matparam_complete(tgt_row,:) = (src_mult + matparam_complete(tgt_row,:) * B_matrix)';
    end

    
    % Fix the incompressibility 1000c1 for the materials without a
    % prescribed K

    % Reading changing materials
    param_names = changing_matrix(4,:);      % Get the parameter names from the 4th row
    materials_number = changing_matrix(1,:);     % Get the values from the 1st row
    
    % Find columns that have 'k' or 'K' (case-insensitive, only exact 'k')
    idx_k = cellfun(@(x) ~isempty(regexpi(x, '^k$', 'once')), param_names);
    
    % Get the corresponding values from the materials numbers
    materials2exclude = cell2mat(materials_number(idx_k));

    % Looping over the elements
    for i = 1:totalNmat

        value2check = i;

        if ismember(value2check, materials2exclude)
            continue
        else
            matparam_complete(i,3) = matparam_complete(i,1)*1000; %1000 times incompressibility
        end
    end

end