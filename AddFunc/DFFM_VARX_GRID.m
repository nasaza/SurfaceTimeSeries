function FPredGrid = DFFM_VARX_GRID(DepV,AddVar1,AddVar2,AddVar3,AddVar4,L_Set,m_Set,h)
% DFFM_VARX_GRID
% Forecasts grid-valued observations by extracting multivariate factors
% directly from station/grid data and forecasting the dependent factor vector
% with a VARX model.
%
% Inputs:
%   DepV    - dependent grid time series, e.g. ozone, N x T
%   AddVar1 - first additional grid time series, N x T
%   AddVar2 - second additional grid time series, N x T
%   AddVar3 - third additional grid time series, N x T
%   AddVar4 - fourth additional grid time series, N x T
%   L_Set   - number of factors: [DepV, AddVar1, AddVar2, AddVar3, AddVar4]
%   m_Set   - number of lags:    [DepV, AddVar1, AddVar2, AddVar3, AddVar4]
%   h       - forecast horizon
%
% Output:
%   FPredGrid - h-step-ahead forecast on the original grid/station locations

%% Step 1: Factor decomposition

PCA_Dep = FactorDecompFTSGrid(DepV,    L_Set(1), 1);
PCA_X1  = FactorDecompFTSGrid(AddVar1, L_Set(2), 1);
PCA_X2  = FactorDecompFTSGrid(AddVar2, L_Set(3), 1);
PCA_X3  = FactorDecompFTSGrid(AddVar3, L_Set(4), 1);
PCA_X4  = FactorDecompFTSGrid(AddVar4, L_Set(5), 1);

%% Step 2: Extract factor scores

FS_Dep = PCA_Dep.scores;
FS_X1  = PCA_X1.scores;
FS_X2  = PCA_X2.scores;
FS_X3  = PCA_X3.scores;
FS_X4  = PCA_X4.scores;

MeanGrid = PCA_Dep.mean;

%% Step 3: Build lagged exogenous regressors for estimation

X1 = lagmatrix(FS_X1,1:m_Set(2));
X2 = lagmatrix(FS_X2,1:m_Set(3));
X3 = lagmatrix(FS_X3,1:m_Set(4));
X4 = lagmatrix(FS_X4,1:m_Set(5));

XX = [X1,X2,X3,X4];

% Remove rows with NaNs created by lagmatrix
valid = all(~isnan(XX),2);

Yest  = FS_Dep(valid,:);
XXest = XX(valid,:);

%% Step 4: Build future exogenous regressors for h-step forecast

XFuture = zeros(h,size(XX,2));

for hh = 1:h

    XF_h = [];

    XF_h = [XF_h, get_future_lags(FS_X1,m_Set(2),hh)];
    XF_h = [XF_h, get_future_lags(FS_X2,m_Set(3),hh)];
    XF_h = [XF_h, get_future_lags(FS_X3,m_Set(4),hh)];
    XF_h = [XF_h, get_future_lags(FS_X4,m_Set(5),hh)];

    XFuture(hh,:) = XF_h;

end

%% Step 5: Estimate and forecast VARX

K     = size(FS_Dep,2);
VARlm = varm(K,m_Set(1));

EsM = estimate(VARlm,Yest,'X',XXest);

Forecast_h = forecast(EsM,h,Yest,'X',XFuture);

%% Step 6: Transform factor forecast back to grid

FPredGrid = PCA_Dep.EigVec(:,1:L_Set(1)) * Forecast_h(h,:)' + MeanGrid;

end


function xLagRow = get_future_lags(FS_X,m,hh)
% Creates the row of lagged exogenous scores needed for forecast step hh.
%
% For h=1, this uses observed lags only.
% For h>1, if a future covariate value is needed but unavailable, the last
% observed covariate score is repeated as a simple persistence convention.

T = size(FS_X,1);

xLagRow = [];

for lag = 1:m

    idx = T + hh - lag;

    if idx <= T
        xLagRow = [xLagRow, FS_X(idx,:)];
    else
        % Future covariate not observed; use last observed value.
        xLagRow = [xLagRow, FS_X(T,:)];
    end

end

end