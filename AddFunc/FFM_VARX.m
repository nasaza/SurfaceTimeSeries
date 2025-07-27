function FuncPredCoeff = FFM_VARX(PCA_Str,PCA_AddVar,L_Set,m_Set,h)

%% Step 1: Create FPCs
    Fscores     = PCA_Str.pcascr;
    coef_eigf   = getcoef(PCA_Str.pcafd);
    coef_mean   = getcoef(PCA_Str.meanfd);
    
%% Step 2: Extarct PCs for additional regressors
    n_add_reg   = length(PCA_AddVar);
    if n_add_reg<1
        disp('No regressors are added. Use FFM_VAR instead')
        return
    end
    
    XX=[];
    for ii =1:n_add_reg  
        XX = [XX, lagmatrix(PCA_AddVar{ii}.pcascr(:,1:L_Set(ii+1)),1:m_Set(ii+1))];
    end

    
%% Step 2: Estimate and forecast L-dimensional VARX(m)     

    L       = L_Set(1);
    Y       = Fscores(:,1:L);
    VARlm   = varm(L,m_Set(1));
    EsM     = estimate(VARlm,Y,'X',XX);

    Forecast_h      = forecast(EsM,h,Y);
    FuncPredCoeff   = (coef_eigf(:,1:L)*(Forecast_h(h,:))')+coef_mean;