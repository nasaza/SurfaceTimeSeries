function FuncPredCoeff   = FRandomForest(PCA_Str,K,h) 

%% 
Fscores     = PCA_Str.pcascr;
coef_eigf   = getcoef(PCA_Str.pcafd);
coef_mean   = getcoef(PCA_Str.meanfd);

%%
[T, ~]      = size(Fscores);   
X_forecast  = zeros(1, K);   

for k = 1:K
    % Create training data:
    % Predictors: all rows from 1 to T-h (using all K columns)
    % Response: the h-step ahead value of series k (rows (1+h) to T)
    X_train = Fscores(1:(T-h), :);
    y_train = Fscores((1+h):T, k);
    
    % Train the random forest regression model using TreeBagger (100 trees)
    rf_model = TreeBagger(100, X_train, y_train, 'Method', 'regression');
    
    % Forecast: use the most recent observation X(T,:) as predictor
    X_forecast(k) = predict(rf_model, Fscores(T, :));
end


FuncPredCoeff   = (coef_eigf(:,1:K)*(X_forecast)')+coef_mean;