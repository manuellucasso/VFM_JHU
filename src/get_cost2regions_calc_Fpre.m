function cost = get_cost2regions_calc_Fpre(path, mymodel, model,edata, x, ...
    p_app,gauss_order,prestress_time,eps,changing_matrix,Normalizer)
    
    global totalRunCount edata_with_Fpre_step ForwardCount
    

    corresponding = changing_matrix(1,:);
    parameters = x .* Normalizer;
    ground_truth_mat = model.matprop;
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


    %Creating material parameters vector to sweep on the PVW calculation
    for ind_mat =1:size(changing_matrix,2)
        col_param = changing_matrix(2,ind_mat);
        if col_param == 2; col_param_prop = 3;else, col_param_prop=col_param; end  
        row_param = changing_matrix(1,ind_mat);
        matparam_sweep(ind_mat,:) = ground_truth_mat(row_param,:);
        matparam_sweep(ind_mat,col_param_prop) = matparam(row_param,col_param); 
    end

    
    % Try-catch block: handle potential failures in FEBio or cost calculation routines.
    try
        % Accumulate prestress calculation for given material parameters updated at every set of parameter.
        % Only update Fpre every 10 evaluations of the cost function
        if ForwardCount==1
            % Always run on the first call, and then every 10th call
            mydir_data = path.data;
            edata = accumulate_Fpre_from_edata(mydir_data, mymodel, gauss_order, prestress_time, matparam, edata,model);
            edata_with_Fpre_step = edata; % Save the new edata
            
            
            % Uncomment for debug:
            fprintf('[Fpre updated at call %d]\n', totalRunCount);
            fprintf('[Number of Forward count call %d]\n', ForwardCount);
            
        else
            % Use the previous Fpre to save on calculation
            edata = edata_with_Fpre_step;
            % Uncomment for debug:
            fprintf('[Fpre reused at call %d]\n', totalRunCount);
        end


        % Compute the cost function using virtual work integrals.
                           [~, ~, ~, cost_func] = calc_virtual_work_variation_integration2(path, mymodel, model, ...
            edata, matparam, matparam_sweep ,p_app, gauss_order,eps,changing_matrix);

        % Check result validity: cost must be finite and non-NaN.
        if isnan(cost_func) || ~isfinite(cost_func)
            error('NaN or Inf encountered in cost function. Triggering restart.');
        end

        % Take logarithm of cost (for scale, stability, or easier optimization).
        cost = log10(cost_func)

        % Debug print (only when totalRunCount == 1 for initialization/testing).
        if totalRunCount == 1
            fprintf('total count %.4f s\n', totalRunCount);
        end

        % Increment global run counter (tracks number of cost calls/runs).
        totalRunCount = totalRunCount+1;
        ForwardCount = ForwardCount+1;

    catch ME
        % Print error message for debugging if something fails in try block.
        fprintf('Restarting due to error: %s\n', ME.message);
        % Rethrow the error so MultiStart/fmincon can handle and possibly restart.
        rethrow(ME);  
    end




    
end