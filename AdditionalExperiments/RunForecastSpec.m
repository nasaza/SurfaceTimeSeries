function Results = RunForecastSpec(d_dyn, L, m, h, saveFolder)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Forecasts Comparison for One Specification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Runs the forecasting comparison for one combination:
%
%   d_dyn = number of lags in cumulative autocovariance operator
%   L     = number of scores/components
%   m     = number of lags in VAR/KNN forecasting step
%   h     = forecast horizon
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Add Libraries

addpath AddFunc
addpath Data

%% Step 1: Read Data
% This step loads seasonally adjusted ozone concentration data observed
% at different stations, the surface data (created in Step 1), and the
% geographic borders of Germany

load('Data\SeasonAdjData'); % Data on the grid seasonaly addjusted
load('Data\FTSs');          % Data given in the surface form (Step 1 has to be executed for that)
ConstrReg = csvread('Data/GeoConstraints/DE_Constraints.csv');

% Projected Coordinates
Threshold = 0.75;

% project data on the plain (Gauss-Krüger Zone 3)
wgs84           = geocrs(4326);
proj            = projcrs(31467);
[LonOz,  LatOz] = projfwd(proj, LatOz, LonOz);
[x, y]          = projfwd(proj, ConstrReg(:,1), ConstrReg(:,2));
ConstrReg       = [y,x];

[XPrecip,YPrecip] = projfwd(proj, YPrecip, XPrecip);
[XWind,  YWind]   = projfwd(proj, YWind, XWind);
[XSun,  YSun]     = projfwd(proj, YSun, XSun);
[XTemp,  YTemp]   = projfwd(proj, YTemp, XTemp);

%% Step 2: Initial values for forecast comparison

H           = 165;  % set forecast horizon
N_Fors      = H+1-h;

[~,T]       = size(OzoneS);
K_max       = 10;
KNN_max     = 50;

q_dyn       = d_dyn;

FEMBasOz    = OzoneBasis;

% Variables to store forecast error surfaces from different predictors

MSurf_MP      = zeros(size(OzoneS,1),N_Fors); % Mean predictor / MF in paper
MSurf_NP      = zeros(size(OzoneS,1),N_Fors); 
MSurf_RF      = zeros(size(OzoneS,1),N_Fors); 
MSurf_FAR1    = zeros(size(OzoneS,1),N_Fors); 

MSurf_PCA_VAR = zeros(size(OzoneS,1),N_Fors);
MSurf_PCA_KNN = zeros(size(OzoneS,1),N_Fors); 

MSurf_DS_VAR  = zeros(size(OzoneS,1),N_Fors); 
MSurf_DS_KNN  = zeros(size(OzoneS,1),N_Fors); 

MSurf_MP_VAR  = zeros(size(OzoneS,1),N_Fors); 
MSurf_MP_KNN  = zeros(size(OzoneS,1),N_Fors); 

% Selected KNN neighbours

KNNs        = zeros(N_Fors,1); % PCA KNN
KNNDS       = zeros(N_Fors,1); % DS KNN
KNNGs       = zeros(N_Fors,1); % MP KNN

% Number of components and lags considered
% IMPORTANT: L and m are the same for static and dynamic scores.

Static_L_set = [L 1 1 1 1]; % Ozone, Sun, Precip, Temp, Wind
Static_m_set = [m 1 1 1 1];

Dyn_L_set    = [L 1 1 1 1];
Dyn_m_set    = [m 1 1 1 1];

%% Step 3: Forecasting with Different Methods

tic

for i = 1:N_Fors

    %% Part of the data to be used for estimation

    nobs_i      = T-H+i-1;
    TrueVal     = OzoneS(:,nobs_i+h);

    %% Static and dynamic scores

    dynstrOz_i  = DynamScoresSurf({OzoneCoef(:,1:nobs_i),OzoneBasis}, K_max, 1, q_dyn); 

    pcastrOz_i  = pca3D({OzoneCoef(:,1:nobs_i),OzoneBasis}, K_max, 1); 

    pcastrS_i   = pca3D({SunCoef(:,1:nobs_i),SunBasis}, K_max, 1); 
    pcastrP_i   = pca3D({PrecipCoef(:,1:nobs_i),PrecipBasis}, K_max, 1); 
    pcastrT_i   = pca3D({TempCoef(:,1:nobs_i),TempBasis}, K_max, 1); 
    pcastrW_i   = pca3D({WindCoef(:,1:nobs_i),WindBasis}, K_max, 1); 

    PCA_addreg  = {pcastrS_i,pcastrP_i,pcastrT_i,pcastrW_i};

    % %% Method 1: Mean Predictor
    % 
    % MSurf_MP(:,i) = TrueVal - mean(OzoneS(:,1:nobs_i),2);
    % 
    % %% Method 2: Naive Predictor
    % 
    % MSurf_NP(:,i) = TrueVal - OzoneS(:,nobs_i);    
    % 
    %% Method 3: FAR(1)

    FuncPredCoeff     = SFAR1(pcastrOz_i,L,h);
    FuncPred          = fd(FuncPredCoeff,FEMBasOz);
    FPointEval        = eval_FEM_fd(LonOz,LatOz,FuncPred);
    MSurf_FAR1(:,i)   = TrueVal - FPointEval;

    %% Method 4: PCA+VAR

    FuncPredCoeff       = FFM_VARX(pcastrOz_i,PCA_addreg,Static_L_set,Static_m_set,h);
    FuncPred            = fd(FuncPredCoeff,FEMBasOz);
    FPointEval          = eval_FEM_fd(LonOz,LatOz,FuncPred);
    MSurf_PCA_VAR(:,i)  = TrueVal - FPointEval;

    % %% Method 5: PCA+KNN
    % %-> Part 1: Estimating number of neighbours
    % 
    % K_NNv        = zeros(KNN_max,1);
    % pcastrKNN_i  = pca3D({OzoneCoef(:,1:nobs_i-h),OzoneBasis}, K_max, 1); 
    % 
    % for knn = 1:KNN_max
    % 
    %     FuncPredCoeff  = FFM_KNN(pcastrKNN_i,Static_L_set(1),Static_m_set(1),knn,h);
    %     FuncPred       = fd(FuncPredCoeff,FEMBasOz);
    %     FPointEval     = eval_FEM_fd(LonOz,LatOz,FuncPred);
    % 
    %     K_NNv(knn,1)   = mean((OzoneS(:,nobs_i)-FPointEval).^2);
    % 
    % end
    % 
    % [~,K_min] = min(K_NNv);
    % KNNs(i,1) = K_min;
    % 
    % %-> Part 2: Final PCA KNN forecast
    % 
    % FuncPredCoeff       = FFM_KNN(pcastrOz_i,Static_L_set(1),Static_m_set(1),K_min,h);
    % FuncPred            = fd(FuncPredCoeff,FEMBasOz);
    % FPointEval          = eval_FEM_fd(LonOz,LatOz,FuncPred);
    % MSurf_PCA_KNN(:,i)  = TrueVal - FPointEval;

    %% Method 6: DS+VAR

    FuncPredCoeff     = FFM_VARX(dynstrOz_i,PCA_addreg,Dyn_L_set,Dyn_m_set,h);
    FuncPred          = fd(FuncPredCoeff,FEMBasOz);
    FPointEval        = eval_FEM_fd(LonOz,LatOz,FuncPred);
    MSurf_DS_VAR(:,i) = TrueVal - FPointEval;

    % %% Method 7: DS+KNN
    % %-> Part 1: Estimating number of neighbours
    % 
    % K_NNv       = zeros(KNN_max,1);
    % 
    % dynstrKNN_i = DynamScoresSurf({OzoneCoef(:,1:nobs_i-h),OzoneBasis}, ...
    %                                K_max, 1, q_dyn);
    % 
    % for knn = 1:KNN_max
    % 
    %     FuncPredCoeff = FFM_KNN(dynstrKNN_i,Dyn_L_set(1),Dyn_m_set(1),knn,h);
    %     FuncPred      = fd(FuncPredCoeff,FEMBasOz);
    %     FPointEval    = eval_FEM_fd(LonOz,LatOz,FuncPred);
    % 
    %     K_NNv(knn,1)  = mean((OzoneS(:,nobs_i)-FPointEval).^2);
    % 
    % end
    % 
    % [~,K_min]  = min(K_NNv);
    % KNNDS(i,1) = K_min;
    % 
    % %-> Part 2: Final DS KNN forecast
    % 
    % FuncPredCoeff     = FFM_KNN(dynstrOz_i,Dyn_L_set(1),Dyn_m_set(1),K_min,h);
    % FuncPred          = fd(FuncPredCoeff,FEMBasOz);
    % FPointEval        = eval_FEM_fd(LonOz,LatOz,FuncPred);
    % MSurf_DS_KNN(:,i) = TrueVal - FPointEval;
    % 
    %% Method 8: MP+VAR

    FPredGrid = DFFM_VARX_GRID(OzoneS(:,1:nobs_i), WindS(:,1:nobs_i), ...
                               SunS(:,1:nobs_i), PrecipS(:,1:nobs_i), ...
                               TempS(:,1:nobs_i), Static_L_set,Static_m_set,h); 

    MSurf_MP_VAR(:,i) = TrueVal - FPredGrid;

    %% Method 9: MP+KNN
    % %-> Part 1: Estimating number of neighbours
    % 
    % K_NNv    = zeros(KNN_max,1);
    % pcastr_i = FactorDecompFTSGrid(OzoneS(:,1:nobs_i-h), K_max, 1);
    % 
    % for knn = 1:KNN_max
    % 
    %     FPredGrid    = DFFM_KNN_GRID(pcastr_i,Static_L_set(1),Static_m_set(1),knn);
    %     K_NNv(knn,1) = mean((OzoneS(:,nobs_i)-FPredGrid).^2);
    % 
    % end
    % 
    % [~,K_min]  = min(K_NNv);
    % KNNGs(i,1) = K_min;
    % 
    % %-> Part 2: Final MP KNN forecast
    % 
    % pcastr_i  = FactorDecompFTSGrid(OzoneS(:,1:nobs_i), K_max, 1);
    % FPredGrid = DFFM_KNN_GRID(pcastr_i,Static_L_set(1),Static_m_set(1),K_min);
    % 
    % MSurf_MP_KNN(:,i) = TrueVal - FPredGrid;

    %% Method 10: DS+Random Forest

    FuncPredCoeff   = FRandomForest(dynstrOz_i,Dyn_L_set(1),h);
    FuncPred        = fd(FuncPredCoeff,FEMBasOz);
    FPointEval      = eval_FEM_fd(LonOz,LatOz,FuncPred);

    MSurf_RF(:,i)   = TrueVal - FPointEval;

    fprintf('Spec d=%d, L=%d, m=%d: forecast %d/%d completed.\n', ...
             d_dyn, L, m, i, N_Fors);

end

RunTime = toc;

%% Analyzing the forecasting performance

MSurf = {MSurf_MP, MSurf_NP, MSurf_FAR1, MSurf_PCA_VAR, MSurf_PCA_KNN, ...
         MSurf_DS_VAR, MSurf_DS_KNN, MSurf_MP_VAR, MSurf_MP_KNN, MSurf_RF};

Models = {'MF','NF','FAR', 'PCA VAR','PCA KNN', 'DS VAR','DS KNN', ...
          'MP VAR','MP KNN', 'DS RF'};

%% Calculate MSE from stored forecast error surfaces

MSE_MP      = mean(MSurf_MP.^2,1)';
MSE_NP      = mean(MSurf_NP.^2,1)';
MSE_FAR1    = mean(MSurf_FAR1.^2,1)';
MSE_PCA_VAR = mean(MSurf_PCA_VAR.^2,1)';
MSE_PCA_KNN = mean(MSurf_PCA_KNN.^2,1)';
MSE_DS_VAR  = mean(MSurf_DS_VAR.^2,1)';
MSE_DS_KNN  = mean(MSurf_DS_KNN.^2,1)';
MSE_MP_VAR  = mean(MSurf_MP_VAR.^2,1)';
MSE_MP_KNN  = mean(MSurf_MP_KNN.^2,1)';
MSE_RF      = mean(MSurf_RF.^2,1)';

MSE_All = [MSE_MP, MSE_NP, MSE_FAR1, MSE_PCA_VAR, MSE_PCA_KNN, MSE_DS_VAR, ...
           MSE_DS_KNN, MSE_MP_VAR, MSE_MP_KNN, MSE_RF];

MeanMSE = mean(MSE_All)

%% Spatial MSE surfaces

MSESurf = cellfun(@(M) mean(M.^2,2), MSurf, 'UniformOutput', false);

%% Save output

Spec.d_dyn = d_dyn;
Spec.L     = L;
Spec.m     = m;
Spec.h     = h;

Results.Spec      = Spec;
Results.Models    = Models;
Results.MSurf     = MSurf;
Results.MSESurf   = MSESurf;
Results.MSE_All   = MSE_All;
Results.MeanMSE   = MeanMSE;
Results.KNNs      = KNNs;
Results.KNNDS     = KNNDS;
Results.KNNGs     = KNNGs;
Results.RunTime   = RunTime;

Results.MSurf_MP      = MSurf_MP;
Results.MSurf_NP      = MSurf_NP;
Results.MSurf_FAR1    = MSurf_FAR1;
Results.MSurf_PCA_VAR = MSurf_PCA_VAR;
Results.MSurf_PCA_KNN = MSurf_PCA_KNN;
Results.MSurf_DS_VAR  = MSurf_DS_VAR;
Results.MSurf_DS_KNN  = MSurf_DS_KNN;
Results.MSurf_MP_VAR  = MSurf_MP_VAR;
Results.MSurf_MP_KNN  = MSurf_MP_KNN;
Results.MSurf_RF      = MSurf_RF;

Results.MSE_MP      = MSE_MP;
Results.MSE_NP      = MSE_NP;
Results.MSE_FAR1    = MSE_FAR1;
Results.MSE_PCA_VAR = MSE_PCA_VAR;
Results.MSE_PCA_KNN = MSE_PCA_KNN;
Results.MSE_DS_VAR  = MSE_DS_VAR;
Results.MSE_DS_KNN  = MSE_DS_KNN;
Results.MSE_MP_VAR  = MSE_MP_VAR;
Results.MSE_MP_KNN  = MSE_MP_KNN;
Results.MSE_RF      = MSE_RF;

if ~exist(saveFolder,'dir')
    mkdir(saveFolder);
end

fileName = sprintf('ForecastSpec_d%d_L%d_m%d_h%d.mat', d_dyn, L, m, h);
save(fullfile(saveFolder,fileName),'Results');

end