function summary = summarize_posterior(posterior_list, nburn, niter, thin, set_zero)
    % Compute posterior summaries (mean, 95% CI, optional thresholding)
    % Inputs:
    %   posterior_list: cell array of p x q matrices
    %   nburn: number of burn-in iterations
    %   niter: number of MCMC iterations after burn-in
    %   thin: thinning interval (default = 1)
    %   set_zero: logical flag to set mean to 0 if CI crosses 0 (default = true)
    % Output:
    %   summary: struct with fields 'mean', 'lower', 'upper'
    
        % Compute indices
        total_iter = nburn + niter;
        idx = nburn+1:thin:total_iter;
        n_keep = numel(idx);
    
        % Get size of one matrix (assume all same size)
        [p, q] = size(posterior_list{idx(1)});
        posterior_array = zeros(p, q, n_keep);
    
        % Stack samples
        for k = 1:n_keep
            posterior_array(:, :, k) = posterior_list{idx(k)};
        end
    
        % Posterior summaries
        mean_est = mean(posterior_array, 3);
        lower_ci = quantile(posterior_array, 0.025, 3);
        upper_ci = quantile(posterior_array, 0.975, 3);
    
        % Set to 0 if CI includes 0
        if set_zero
            zero_mask = (lower_ci .* upper_ci < 0);
            mean_est(zero_mask) = 0;
        end
    
        % Return result as struct
        summary.mean = mean_est;
        summary.lower = lower_ci;
        summary.upper = upper_ci;
    end
    