function FPredGrid = DFFM_VARX_GRID(DepV,AddVar1,AddVar2,AddVar3,AddVar4,L_Set,m_Set,h)
% This procedure performs forecasting based on factros extracted from 
%  several time series; Based on Aue et al (2015) Section 4

%% Step 1: Factro Decomposition

    PCA_Dep     = FactorDecompFTSGrid(DepV, L_Set(1), 1);
    PCA_X1      = FactorDecompFTSGrid(AddVar1, L_Set(2), 1);
    PCA_X2      = FactorDecompFTSGrid(AddVar2, L_Set(3), 1);
    PCA_X3      = FactorDecompFTSGrid(AddVar3, L_Set(4), 1);
    PCA_X4      = FactorDecompFTSGrid(AddVar4, L_Set(5), 1);

%% Step 1: Create FPCs
    FS_Dep      = PCA_Dep.scores;
    FS_X1       = PCA_X1.scores;
    FS_X2       = PCA_X2.scores;
    FS_X3       = PCA_X3.scores;
    FS_X4       = PCA_X4.scores;
    
    X1          = lagmatrix(FS_X1,1:m_Set(2));
    X2          = lagmatrix(FS_X2,1:m_Set(3));
    X3          = lagmatrix(FS_X3,1:m_Set(4));
    X4          = lagmatrix(FS_X4,1:m_Set(5));
    
    XX          = [X1,X2,X3,X4];
    
    MeanGrid    = PCA_Dep.mean;
    
    
%% Step 2: Estimate and forecast L-dimensional VAR(m) 
    K           = size(FS_Dep,2);
    VARlm       = varm(K,max(m_Set));
    EsM         = estimate(VARlm,FS_Dep,'X',XX);

    Forecats_h  = forecast(EsM,h,FS_Dep);
    FPredGrid   = (PCA_Dep.EigVec(:,1:L_Set(1))*Forecats_h(h,:)') + MeanGrid;