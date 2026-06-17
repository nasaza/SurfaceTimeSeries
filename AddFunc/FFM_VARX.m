function FuncPredCoeff = FFM_VARX(PCA_Str,PCA_AddVar,L_Set,m_Set,h)
% FFM_VARX
% Forecasts functional principal component scores using VAR or VARX.
%
% If PCA_AddVar is empty, the function estimates a pure VAR model for the
% response scores. If PCA_AddVar is non-empty, additional lagged score
% regressors are included as exogenous variables.
%
% Inputs:
%   PCA_Str    - structure with fields pcascr, pcafd, meanfd
%   PCA_AddVar - cell array of additional PCA structures
%   L_Set      - number of components [response, regressor 1, ..., regressor k]
%   m_Set      - lag orders [response, regressor 1, ..., regressor k]
%   h          - forecast horizon
%
% Output:
%   FuncPredCoeff - coefficient vector of the h-step-ahead forecast surface

%% Step 1: Create response score matrix

Fscores   = PCA_Str.pcascr;
coef_eigf = getcoef(PCA_Str.pcafd);
coef_mean = getcoef(PCA_Str.meanfd);

L = L_Set(1);
Y = Fscores(:,1:L);

%% Step 2: Pure VAR case

VARlm = varm(L,m_Set(1));

if nargin < 2 || isempty(PCA_AddVar)

    EsM        = estimate(VARlm,Y);
    Forecast_h = forecast(EsM,h,Y);

else

    %% Step 3: Extract PCs for additional regressors

    n_add_reg = length(PCA_AddVar);

    XX = [];

    for ii = 1:n_add_reg

        Xi_scores = PCA_AddVar{ii}.pcascr(:,1:L_Set(ii+1));

        XX = [XX, lagmatrix(Xi_scores,1:m_Set(ii+1))];

    end

    %% Step 4: Remove rows with missing lagged regressors

    valid = all(~isnan(XX),2);

    Yest  = Y(valid,:);
    XXest = XX(valid,:);

    %% Step 5: Construct future exogenous regressors correctly

    XFuture = zeros(h,size(XX,2));

    for hh = 1:h

        XF_h = [];

        for ii = 1:n_add_reg

            Xi_scores = PCA_AddVar{ii}.pcascr(:,1:L_Set(ii+1));

            XF_h = [XF_h, get_future_lags_scores(Xi_scores,m_Set(ii+1),hh)];

        end

        XFuture(hh,:) = XF_h;

    end

    %% Step 6: Estimate and forecast VARX

    EsM        = estimate(VARlm,Yest,'X',XXest);
    Forecast_h = forecast(EsM,h,Yest,'X',XFuture);

end

%% Step 7: Transform forecasted scores back to functional coefficients

FuncPredCoeff = coef_eigf(:,1:L) * Forecast_h(h,:)' + coef_mean;

end


function xLagRow = get_future_lags_scores(Xscores,m,hh)
% get_future_lags_scores
% Creates the future row of lagged exogenous score regressors.
%
% In estimation, lagmatrix(Xscores,1:m) gives regressors
% X_{t-1}, X_{t-2}, ..., X_{t-m} for equation Y_t.
%
% Therefore, for forecasting Y_{T+hh}, we need
% X_{T+hh-1}, X_{T+hh-2}, ..., X_{T+hh-m}.
%
% If a required future covariate value is unavailable, we use the last
% observed covariate score as a simple persistence convention.

T = size(Xscores,1);

xLagRow = [];

for lag = 1:m

    idx = T + hh - lag;

    if idx <= T
        xLagRow = [xLagRow, Xscores(idx,:)];
    else
        % Future covariate score is not observed.
        % Use last observed value as a persistence approximation.
        xLagRow = [xLagRow, Xscores(T,:)];
    end

end

end