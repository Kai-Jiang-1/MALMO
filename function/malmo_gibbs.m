% MATLAB Version of MALMO Gibbs Sampler
% Function definition
function result = malmo_gibbs(Y, X, nburn, niter, seed_num)
    % Gibbs Sampling Function for MALMO
    % Inputs:
    %   Y: n x q response matrix (log metabolite values); missing entries must be coded as 0
    %   X: n x p predictor matrix (microbiome abundances, centered and scaled)
    %   nburn: Number of burn in iterations
    %   niter: Number of iterations kept after burn in
    %   seed_num: Seed number for reproducibility

    % Set random seed
    rng(seed_num);

    % Get dimensions
    [n, q] = size(Y);
    p = size(X, 2);

    % Initialize result storage
    B_result = cell(nburn + niter, 1);
    A_result = cell(nburn + niter, 1);
    C_result = cell(nburn + niter, 1);
    Sigma_result = cell(nburn + niter, 1);
    beta0_result = cell(nburn + niter, 1);
    lambda2_result = cell(nburn + niter, 1);
    tau2_result = cell(nburn + niter, 1);
    psi_result = cell(nburn + niter, 1);
    nu_result = cell(nburn + niter, 1);
    Z_result = cell(nburn + niter, 1);
    Y_mvi_result = cell(nburn + niter, 1);
    theta_result = zeros(nburn + niter, 1);


    % ------------------- Initialize Parameters -------------------
    % nu: Inverse Gamma (shape = 0.5, scale = 1)
    nu = 1 ./ gamrnd(0.5, 1, [p, q]); % Inverse Gamma in MATLAB: 1/gamrnd(a, 1/scale)
    nu_result{1} = nu;

    % lambda2: Inverse Gamma (shape = 0.5, scale = 1/nu_jh)
    lambda2 = 1 ./ gamrnd(0.5, nu, [p, q]);
    lambda2_result{1} = lambda2;

    % psi: Inverse Gamma (shape = 0.5, scale = 1)
    psi = 1 ./ gamrnd(0.5, 1, [q, 1]);
    psi_result{1} = psi;

    % tau2: Inverse Gamma (shape = 0.5, scale = 1/psi)
    tau2 = 1 ./ gamrnd(0.5, psi, [q, 1]);
    tau2_result{1} = tau2;

    % Sigma: Diagonal matrix with sigma_h^2 = 1 initially
    sigma2_values = ones(q,1);
    Sigma = diag(sigma2_values);
    Sigma_result{1} = Sigma;

    % A: a_hk ~ N(0, 1), h, k = 1,...,q
    A = randn(q, q); % Normal distribution with mean 0, std 1
    A_result{1} = A;

    % beta0: beta0_h ~ N(0, Sigma_hh)
    beta0 = randn(q, 1) .* sqrt(diag(Sigma));
    beta0_result{1} = beta0;

    % Missingness indicator (R: 0 = missing, 1 = observed)
    % Note: any entry exactly equal to 0 is treated as missing
    R = zeros(n, q);
    R(Y ~= 0) = 1;

    % Limit of Detection (xi): smallest observed value in Y, shared by all metabolites
    LOD = min(Y(Y ~= 0));
    % theta: Uniform(0, 1)
    theta = unifrnd(0, 1);
    theta_result(1) = theta;


    % ------------------- Gibbs Sampling Loop -------------------
    for i = 1:(nburn + niter)
        % Sample B (exact fast sampler of Bhattacharya, Chakraborty and Mallick, 2016)
        % Y_tilde = Y - 1_n * beta0'
        Y_tilde = Y - ones(n, 1) * beta0';
        y = reshape(Y_tilde', [], 1); % Convert Y_tilde to a column vector

        % Sigma_tilde = I_n kron Sigma
        Sigma_diag = diag(Sigma);
        inv_sqrt_Sigma_diag = 1 ./ sqrt(Sigma_diag);
        inv_sqrt_diag_vec = repmat(inv_sqrt_Sigma_diag, n, 1);
        Sigma_tilde_inv_sqrt = diag(inv_sqrt_diag_vec);

        % Then use:
        % X_tilde = Sigma_tilde_inv_sqrt * kron(X, A);
        X_tilde = Sigma_tilde_inv_sqrt * kron(X, A);

        % y_tilde = (Sigma_tilde)^(-1/2) * y
        y_tilde = Sigma_tilde_inv_sqrt * y;

        % Create Lambda, a diagonal matrix of variances
        Lambda_vec = reshape((lambda2 .* tau2')', [], 1);
        Lambda = diag(Lambda_vec);

        % Step 1: Sample u and delta
        u = randn(p * q, 1) .* sqrt(diag(Lambda)); % Same as sqrt(Lambda_vec)
        delta = randn(n * q, 1);

        % Set v = X_tilde * u + delta
        v = X_tilde * u + delta;

        % Solve for w using Cholesky decomposition for stability and speed
        middle_mat = X_tilde * Lambda * X_tilde' + eye(n * q);
        % middle_mat_chol = chol(middle_mat, 'upper');
        % w = middle_mat_chol \ (middle_mat_chol' \ (y_tilde - v));

        try
            middle_mat_chol = chol(middle_mat, 'upper');
            w = middle_mat_chol \ (middle_mat_chol' \ (y_tilde - v));
        catch
            warning_state = warning('off', 'MATLAB:nearlySingularMatrix');
            w = linsolve(middle_mat, (y_tilde - v));
            warning(warning_state);
        end

        % Set beta = u + Lambda * X_tilde' * w
        beta = u + Lambda * X_tilde' * w;

        % Reshape beta into a matrix B (p x q)
        % The transpose here correctly handles R's `byrow = TRUE`
        B = reshape(beta, q, p)';
        B_result{i} = B;

        %% Sample matrix A
        X_star = Sigma_tilde_inv_sqrt * kron(X * B, eye(q));
        omega_a = X_star' * X_star + eye(q^2);
        omega_a_chol = chol(omega_a, 'lower');
        v_a = omega_a_chol \ (X_star' * y_tilde);
        m_a = omega_a_chol' \ v_a;
        z_a = randn(q^2, 1);
        w_a = omega_a_chol' \ z_a;
        a = m_a + w_a;
        A = reshape(a, q, q);
        A_result{i} = A;

        % Calculate and save the final coefficient matrix C = B * A'
        C = B * A';
        C_result{i} = C;

        % Sample local shrinkage parameters
        % lambda2_jh ~ IG(1, (B^2_jh / (2 * tau2_h)) + (1 / nu_jh))
        scale_for_lambda2 = (B.^2) ./ (2 * repmat(tau2', p, 1)) + (1 ./ nu);
        lambda2 = 1 ./ gamrnd(1, 1 ./ scale_for_lambda2, p, q);
        lambda2_result{i} = lambda2;

        % nu_jh ~ IG(1, 1 + 1 / lambda2_jh)
        scale_for_nu = 1 + (1 ./ lambda2);
        nu = 1 ./ gamrnd(1, 1 ./ scale_for_nu, p, q);
        nu_result{i} = nu;

        % Sample global shrinkage parameters
        % tau2_h ~ IG((p + 1)/2, 0.5 * sum_j (B_jh^2 / lambda2_jh) + (1 / psi_h))
        scale_for_tau2 = 0.5 * sum(B.^2 ./ lambda2, 1)' + (1 ./ psi);
        tau2 = 1 ./ gamrnd((p + 1) / 2, 1 ./ scale_for_tau2, q, 1);
        tau2_result{i} = tau2;

        % psi_h ~ IG(1, 1 + 1 / tau2_h)
        scale_for_psi = 1 + (1 ./ tau2);
        psi = 1 ./ gamrnd(1, 1 ./ scale_for_psi, q, 1);
        psi_result{i} = psi;

        % Sample beta0 and sigma
        % Step 1: Residual
        Y_star = Y - X * C;

        % Step 2: Mean of beta0
        sum_Y_star = sum(Y_star, 1);
        mean_beta0 = sum_Y_star' / (n + 1); % mean_beta0 will be q x 1

        % Step 3: Cov of beta0
        cov_beta0 = diag(Sigma) ./ (n + 1);
        beta0 = mean_beta0 + sqrt(cov_beta0) .* randn(q, 1);
        beta0_result{i} = beta0;

        % Sigma
        % sigma2_h ~ IG((n + 1)/2, (SSR_h + beta0_h^2)/2), since beta0_h ~ N(0, sigma2_h)
        Residuals = Y - ones(n, 1) * beta0' - X * C;
        shape_sigma2 = (n + 1) / 2;
        scale_sigma2 = 0.5 * (sum(Residuals.^2, 1)' + beta0.^2);
        sigma2_values = 1 ./ gamrnd(shape_sigma2 * ones(q, 1), 1 ./ scale_sigma2);
        Sigma = diag(sigma2_values);
        Sigma_result{i} = Sigma;

        mu_matrix = ones(n, 1) * beta0' + X * C;
        sigma_sds = repmat(sqrt(diag(Sigma))', n, 1);

        % Identify missing values (R == 0)
        is_miss = (R == 0); % is_miss will be a logical n x q matrix

        phi = nan(n, q);
        phi(is_miss) = normcdf(LOD, mu_matrix(is_miss), sigma_sds(is_miss));

        phi_value = nan(n, q);
        phi_value(is_miss) = 1 ./ (1 + theta * ((1 - phi(is_miss)) ./ phi(is_miss)));

        % Step 3: Sample Z
        % z = 2 observed
        % z = 1 missing and lower than LOD
        % z = 0 missing and higher than LOD
        Z = 2 * ones(n, q); % Initialize all as observed (2)
        Z(is_miss) = binornd(1, phi_value(is_miss));
        Z_result{i} = Z;

        Y_mvi = nan(n, q);

        mask_low = (Z == 1); % less than lod
        if any(mask_low, 'all')
            Y_mvi(mask_low) = truncnormrnd(-Inf, LOD, mu_matrix(mask_low), sigma_sds(mask_low));
        end

        mask_high = (Z == 0); % greater than lod
        if any(mask_high, 'all')
            Y_mvi(mask_high) = truncnormrnd(LOD, Inf, mu_matrix(mask_high), sigma_sds(mask_high));
        end

        Y_mvi(Z == 2) = Y(Z == 2);
        Y_mvi_result{i} = Y_mvi;

        theta = betarnd(sum(Z == 0, 'all') + 1, sum(Z == 2, 'all') + 1);
        theta_result(i) = theta;
        Y = Y_mvi;

        % --- Display iteration progress ---
        if mod(i, 1000) == 0
            fprintf('Iteration: %d\n', i);
        end
    end

    % Always return these results
    result.B_result = B_result;
    result.A_result = A_result;
    result.C_result = C_result;
    result.Sigma_result = Sigma_result;
    result.beta0_result = beta0_result;
    result.lambda2_result = lambda2_result;
    result.tau2_result = tau2_result;
    result.psi_result = psi_result;
    result.nu_result = nu_result;
    result.Z_result = Z_result;
    result.Y_mvi_result = Y_mvi_result;
    result.theta_result = theta_result;
end