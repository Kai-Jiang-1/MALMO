
# MALMO

Multivariate Adaptive Linear Model with Missing Outcomes

**Author**: Kai Jiang, Christine B. Peterson

## Overview

This repository provides the implementation of the **Multivariate Adaptive Linear Model with Missing Outcomes (MALMO)**.

The proposed framework enables:

- **Joint modeling** of high-dimensional multivariate outcomes
- **Variable selection** via horseshoe prior
- **Missing value imputation** for multivariate outcomes
- **Bayesian inference** via Gibbs sampling

## Repository Structure

- `simulation/` – Scripts for generating simulated datasets
- `function/` – Contains MALMO model and related functions

## Example

### Step 1 - Load Simulated Data
```matlab
filename = sprintf('simulation/%s_n%d_q%d_p%d_snr%d_prop%g.mat', ...
                   'mnar0.33', 50, 3, 100, 1, 0.2);
load(filename);
whos
```

### Step 2 - Scale and Split Data
```matlab
mu_X = mean(X, 1);
sd_X = std(X, 0, 1);
X = (X - mu_X) ./ sd_X;

rng(1996);
n_train = floor(0.7 * size(X,1));
perm_idx = randperm(size(X,1));
train_idx = sort(perm_idx(1:n_train));
test_idx  = sort(perm_idx(n_train+1:end));

X_train = X(train_idx, :);
X_test  = X(test_idx, :);
Y_train = Y_miss(train_idx, :);
Y_test  = Y_miss(test_idx, :);
```

### Step 3 - Run MALMO
```matlab
nburn = 5000;
niter = 10000;
seed_num = 123;
result = malmo_gibbs(Y_train, X_train, nburn, niter, seed_num);
```

### Step 4 - Post Processing
```matlab
thin = 10;

% Coefficients
C_summary = summarize_posterior(result.C_result, nburn, niter, thin, true);
C_est = C_summary.mean;

% Y_mvi
Y_mvi_summary = summarize_posterior(result.Y_mvi_result, nburn, niter, thin, false);
```

### Step 5 - Evaluate Results
```matlab
C_select_row = double(any(C_est ~= 0, 2));
C_idx_row    = double(any(C_idx ~= 0, 2));
TP_row  = sum(C_select_row == 1 & C_idx_row == 1);
FP_row  = sum(C_select_row == 1 & C_idx_row == 0);
FN_row  = sum(C_select_row == 0 & C_idx_row == 1);
TN_row  = sum(C_select_row == 0 & C_idx_row == 0);
TPR_row = TP_row / (TP_row + FN_row);
FPR_row = FP_row / (FP_row + TN_row);
```

