function FuncPredCoeff = FFM_KNN(PCA_Str,L,m,KNN,h)

%% Step 1: PCA and initial values

    Fscores     = PCA_Str.pcascr;
    coef_eigf   = getcoef(PCA_Str.pcafd);
    coef_mean   = getcoef(PCA_Str.meanfd);
    FeatureSc   = Fscores(end-m+1:end,1:L); % define feature scores: "last" observations including lags
    LearningSc  = Fscores(1:end-h,1:L);     % define lerning sample scores
    [T,~]       = size(LearningSc);         % Overall number of available observations
    SortingNN   = zeros(T+1-m,2);           % Matrix to collect distance to NN
    weights     = PCA_Str.varprop;          % create for scores according to their eigenvalue
    weights     = weights(1:L);
    WW          = repmat(weights',m,1);
    
%% Step 2: Calculate distances to NN

    for ii=m:T
        SortingNN(ii-m+1,1) = ii;
        SortingNN(ii-m+1,2) = sqrt(1/m*sum(sum((WW.*(FeatureSc - LearningSc(ii-m+1:ii,:))).^2)));
    end
    
 %% Step 3: Selectin NN and constracting a forecast
    
    dist_sort       = sortrows(SortingNN,2);        % sort distances by increasing value
    idx             = dist_sort(1:KNN,1);           % retain indices of nearest neighbors
    dist_s          = dist_sort(1:KNN,2)+0.000001;  % retain distance to nearest neighbors
    score_forecast  = zeros(1,L);                   % initialize SCORE_FORECAST
    score_mat       = Fscores(:,1:L);
    sum_dist        = sum(1./dist_s);
    for ii=1:KNN
        score_forecast = score_forecast + ((1/dist_s(ii))/ sum_dist .* score_mat(idx(ii)+1,:));
    end
    FuncPredCoeff   = (score_forecast * coef_eigf(:,1:L)')'+coef_mean;
    
