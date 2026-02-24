function [matparam,matparam_sweep] = SweepMatrix(model,changing_matrix,x,Normalizer)
    % Compute Cauchy stress for material in 'model' at index 'mat_idx'
    % F: 3x3 deformation gradient
    % model: struct from your parser (see above)
    % mat_idx: index of material to use
           
    ground_truth_mat = model.matprop;
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
        
        switch mat_type
            case "neo-hookean"
                switch lower(mat_change)
                    case 'mu'
                        prop=1;
                    case 'k'
                        prop=2;
                end
            case "mooney-rivlin"
                switch lower(mat_change)
                    case 'c1'
                        prop=1;
                    case 'c2'
                        prop=2;
                    case 'k'
                        prop=3;
                end
            case "ti-mooney-rivlin"
                switch lower(mat_change)
                    case 'c1'
                        prop=1;
                    case 'c2'
                        prop=2;
                    case 'k'
                        prop=3;
                    case 'c3'
                        prop=4;
                    case 'c4'
                        prop=5;
                    case 'c5'
                        prop=6;
                    case 'flam'
                        prop=7;
                end
            case 'hgo-unconstrained'
                switch lower(mat_change)
                    case 'c'
                        prop=1;
                    case 'k1'
                        prop=2;
                    case 'k2'
                        prop=4;
                    case 'k'
                        prop=3;
                    case 'kappa'
                        prop=5;
                    case 'gamma'
                        prop=6;
                end
        end

    
    % Changing the value in the material parameter sweep matrix
        col_param_prop = prop;
        row_param = changing_matrix{1,col}; 
        col_param = changing_matrix{2,col};  
        matparam_sweep(col,:) = ground_truth_mat(row_param,:);
        matparam_sweep(col,col_param_prop) = matparam(row_param,col_param);
    
    end

end