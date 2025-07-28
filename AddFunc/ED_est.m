function L = ED_est(Eig_values)
% Calculates an estimator of the number of factors following Onatski (2010)

%% Initial values
K_max = length(Eig_values);
y     = Eig_values(K_max-4:K_max);
X     = (K_max-1:K_max+3)';
X     = X.^(2/3);
d     = 2*abs((X'*X)^(-1)*(X'*y));

%% Compute Estimate of L
v_diff = Eig_values(1:K_max-1)-Eig_values(2:K_max);
L      = max(find(v_diff>=d));