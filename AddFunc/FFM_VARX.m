function FuncPredCoeff = FFM_VARX(PCA_Str,PCA_AddVar,L_Set,m_Set,h)

%% Step 1: Create FPCs
    Fscores     = PCA_Str.pcascr;
    coef_eigf   = getcoef(PCA_Str.pcafd);
    coef_mean   = getcoef(PCA_Str.meanfd);
    L           = L_Set(1);
    Y           = Fscores(:,1:L);
    
%% Step 2: Extarct PCs for additional regressors
    XX = [];
    if nargin >= 2 && ~isempty(PCA_AddVar)
        n_add_reg = length(PCA_AddVar);
        for ii = 1:n_add_reg
            XX = [XX, ...
                  lagmatrix(PCA_AddVar{ii}.pcascr(:,1:L_Set(ii+1)), ...
                            1:m_Set(ii+1))];
        end
    end
    
%% Step 2: Estimate and forecast L-dimensional VARX(m)     

    VARlm   = varm(L,m_Set(1));
    if isempty(XX)   % Pure VAR case
        EsM         = estimate(VARlm,Y);
        Forecast_h  = forecast(EsM,h,Y);
    else % VARX case
        valid       = all(~isnan(XX),2);    
        Yest        = Y(valid,:);
        XXest       = XX(valid,:);
        EsM         = estimate(VARlm,Yest,'X',XXest);
        XFuture     = repmat(XXest(end,:),h,1);
        Forecast_h  = forecast(EsM,h,Yest,'X',XFuture);
    end

    FuncPredCoeff   = (coef_eigf(:,1:L)*(Forecast_h(h,:))')+coef_mean;
    
end