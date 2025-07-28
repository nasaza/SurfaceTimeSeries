function [OS,ANH] = InformCriteria(PCA_Str,K_max,p_max)
% produses estimates from two information criteria Aue et al. (2015) (in
% ANH) and Otto and Salish(2021) (in OS).
%% Step 1: Create FPCs
Fscores         = PCA_Str.pcascr;
Feigval         = PCA_Str.values;
[T,~]           = size(Fscores);

%% Step 2: Calculate MSM criteria for each K and p
Lm_matANH          = zeros(K_max,p_max);
MSEf               = zeros(K_max,p_max);
Lm_matOSBIC        = zeros(K_max,p_max);
Lm_matOSHQ         = zeros(K_max,p_max);
for l=1:K_max
    for m =1:p_max
        VARlm           = varm(l,m);
        [EsM]           = estimate(VARlm,Fscores(:,1:l));
        Lm_matANH(l,m)  = (T+l*m)/(T-l*m)*trace(EsM.Covariance)+sum(Feigval(l+1:end));
        MSEf(l,m)       = trace(EsM.Covariance)+sum(Feigval(l+1:end));
        Lm_matOSBIC(l,m)= log(trace(EsM.Covariance)+sum(Feigval(l+1:end)))+l*m*log(T)/T;
        Lm_matOSHQ(l,m) = log(trace(EsM.Covariance)+sum(Feigval(l+1:end)))+2*l*m*log(log(T))/T;
    end    
end

%% Slecting estimates of K and p
[Min_val,L_vect]= min(Lm_matANH);
[~,p_est]       = min(Min_val);
L_est           = L_vect(p_est);
ANH.L           = L_est;
ANH.m           = p_est;

[Min_val,L_vect]= min(Lm_matOSBIC);
[~,p_est]       = min(Min_val);
L_est           = L_vect(p_est);
OS.L            = L_est;
OS.m            = p_est;

[Min_val,L_vect]= min(Lm_matOSHQ);
[~,p_est]       = min(Min_val);
L_est           = L_vect(p_est);
OS.HQ_L         = L_est;
OS.HQ_m         = p_est;

OS.MSE          = MSEf;
