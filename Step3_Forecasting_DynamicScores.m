%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Forecasts Comparison
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
clc;

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
% This step sets up initial values for the forecast comparison: the number of
% observations saved and used for forecast comparison (H), the forecasting
% horizon (h), the maximum number of lags considered (p_max), and the maximum
% number of neighbors in the KNN method (KNN_max). Additionally, it creates "empty"
% variables where forecasts will be stored (e.g., MSE_... and Msurf...).

H           = 165;  %set forecast horizon
h           = 1;    % h-step ahead
N_Fors      = H+1-h;

[~,T]       = size(OzoneS);
K_max       = 10;
p_max       = 5;   % Maximum nuber of lag considered in the analysis
KNN_max     = 50;  % Maximum number of KNN in the calibraiton
KNNs        = zeros(H,2);
q_dyn       = 3;

FEMBasOz    = OzoneBasis;

% Variables to store MSE surfaces from different predictors
MSurf_MP      = zeros(size(OzoneS,1),(H+1-h)); 
MSurf_NP      = zeros(size(OzoneS,1),(H+1-h)); 
MSurf_RF      = zeros(size(OzoneS,1),(H+1-h)); 
MSurf_FAR1    = zeros(size(OzoneS,1),(H+1-h)); 
MSurf_PCA_VAR = zeros(size(OzoneS,1),(H+1-h));
MSurf_PCA_KNN = zeros(size(OzoneS,1),(H+1-h)); 
MSurf_DS_VAR  = zeros(size(OzoneS,1),(H+1-h)); 
MSurf_DS_KNN  = zeros(size(OzoneS,1),(H+1-h)); 
MSurf_MP_VAR  = zeros(size(OzoneS,1),(H+1-h)); 
MSurf_MP_KNN  = zeros(size(OzoneS,1),(H+1-h)); 

% Number of components and lags considered 
Static_L_set = [3 1 1 1 1]; % Ozone, Sun, Precip, Temp, Wind
Static_m_set = [2 1 1 1 1];
Dyn_L_set    = [4      1 1 1 1]; 
Dyn_m_set    = [3 1 1 1 1];


%% Step 3: Forecasting with Different Methods
% Calculates forecasts for horizon h over N_Fors iterations, updating the sample each time
% and saving the forecast outputs in the corresponding variables

tic

for i=1:(H+1-h)    
%% Part of the data to be used for estimation
    nobs_i      = T-H+i-1;
    TrueVal     = OzoneS(:,nobs_i+h);    
    OzoneFTSi   = fd(OzoneCoef(:,1:nobs_i),OzoneBasis);
    
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
%     MSurf_MP(:,i) = TrueVal-mean(OzoneS(:,1:nobs_i),2);
% 
% %% Method 2: Naive Predictor
% 
%     MSurf_NP(:,i) = TrueVal-OzoneS(:,nobs_i);    
% 
% %% Method 3: FAR(1)
% 
%     FuncPredCoeff   = SFAR1(pcastrOz_i,Static_L(1),h);
%     FuncPred        = fd(FuncPredCoeff,FEMBasOz);
%     FPointEval      = eval_FEM_fd(LonOz,LatOz,FuncPred);
%     MSurf_FAR1(:,i)     = TrueVal-FPointEval;
% 
%% Method 4: PCA+VAR

    FuncPredCoeff   = FFM_VARX(pcastrOz_i,PCA_addreg,Static_L_set,Static_m_set,h);
    FuncPred        = fd(FuncPredCoeff,FEMBasOz);
    FPointEval      = eval_FEM_fd(LonOz,LatOz,FuncPred);
    MSurf_PCA_VAR(:,i)     = TrueVal-FPointEval;

% %% Method 5: PCA+KNN
% %-> Part 1: Estimating number of neigbours
% 
%     K_NNv        = zeros(KNN_max,1) ;
%     pcastrKNN_i  = pca3D({OzoneCoef(:,1:nobs_i-h),OzoneBasis}, K_max, 1); 
%     for knn = 1:KNN_max
%         FuncPredCoeff  = FFM_KNN(pcastrKNN_i,Static_L_set(1),Static_m_set(1),knn,h);
%         FuncPred       = fd(FuncPredCoeff,FEMBasOz);
%         FPointEval     = eval_FEM_fd(LonOz,LatOz,FuncPred);
%         K_NNv(knn,1)   = mean((OzoneS(:,nobs_i)-FPointEval).^2);
%     end
%     [~,K_min]     = min(K_NNv);
%     KNNs(i,:)     = K_min;
% %-> Part 2    
%     FuncPredCoeff  = FFM_KNN(pcastrOz_i,Static_L_set(1),Static_m_set(1),K_min,h);
%     FuncPred       = fd(FuncPredCoeff,FEMBasOz);
%     FPointEval     = eval_FEM_fd(LonOz,LatOz,FuncPred);
%     MSurf_PCA_KNN(:,i)    = TrueVal-FPointEval;

%% Method 6: DS+VAR

    FuncPredCoeff = FFM_VARX(dynstrOz_i,PCA_addreg,Dyn_L_set,Dyn_m_set,h);
    FuncPred      = fd(FuncPredCoeff,FEMBasOz);
    FPointEval    = eval_FEM_fd(LonOz,LatOz,FuncPred);
    MSurf_DS_VAR(:,i) = TrueVal - FPointEval;

% %% Method 7: DS+KNN
%     %-> Part 1: Estimating number of neighbours
% 
%     K_NNv       = zeros(KNN_max,1);
%     dynstrKNN_i = DynamScoresSurf({OzoneCoef(:,1:nobs_i-h),OzoneBasis}, ...
%                                    K_max, 1, q_dyn);
%     for knn = 1:KNN_max
%         FuncPredCoeff = FFM_KNN(dynstrKNN_i,Dyn_L_set(1),Dyn_m_set(1),knn,h);
%         FuncPred      = fd(FuncPredCoeff,FEMBasOz);
%         FPointEval    = eval_FEM_fd(LonOz,LatOz,FuncPred);
%         K_NNv(knn,1)  = mean((OzoneS(:,nobs_i)-FPointEval).^2);
% 
%     end
%     [~,K_min]   = min(K_NNv);
%     KNNDS(i,1)  = K_min;
% 
%     %-> Part 2: Final DS KNN forecast
%     FuncPredCoeff = FFM_KNN(dynstrOz_i,Dyn_L_set(1),Dyn_m_set(1),K_min,h);
%     FuncPred      = fd(FuncPredCoeff,FEMBasOz);
%     FPointEval    = eval_FEM_fd(LonOz,LatOz,FuncPred);
%     MSurf_DS_KNN(:,i) = TrueVal - FPointEval;

    
 % %% Method 8: MP+VAR
 %    FPredGrid = DFFM_VARX_GRID(OzoneS(:,1:nobs_i), WindS(:,1:nobs_i), ...
 %                               SunS(:,1:nobs_i), PrecipS(:,1:nobs_i), ...
 %                               TempS(:,1:nobs_i), Static_L_set,Static_m_set,h); 
 %    MSurf_MP_VAR(:,i) = TrueVal - FPredGrid;
 % 
 %    %% Method 9: MP+KNN
 %    %-> Part 1: Estimating number of neighbours
 %    K_NNv       = zeros(KNN_max,1);
 %    pcastr_i    = FactorDecompFTSGrid(OzoneS(:,1:nobs_i-h), K_max, 1);
 % 
 %    for knn = 1:KNN_max
 %        FPredGrid    = DFFM_KNN_GRID(pcastr_i,Static_L_set(1),Static_m_set(1),knn);
 %        K_NNv(knn,1) = mean((OzoneS(:,nobs_i)-FPredGrid).^2);
 %    end
 % 
 %    [~,K_min]  = min(K_NNv);
 %    KNNGs(i,1) = K_min;
 % 
 %    %-> Part 2: Final MP KNN forecast
 % 
 %    pcastr_i  = FactorDecompFTSGrid(OzoneS(:,1:nobs_i), K_max, 1);
 %    FPredGrid = DFFM_KNN_GRID(pcastr_i,Static_L_set(1),Static_m_set(1),K_min);
 %    MSurf_MP_KNN(:,i) = TrueVal - FPredGrid;
 % 
 % 

%% Method 10: DS+Random Forest

    % FuncPredCoeff   = FRandomForest(dynstrOz_i,Dyn_L_set(1),h);
    % FuncPred        = fd(FuncPredCoeff,FEMBasOz);
    % FPointEval      = eval_FEM_fd(LonOz,LatOz,FuncPred);
    % MSurf_RF(:,i)   = TrueVal - FPointEval;


    i
end
toc

%% Analyzing the forecasting performance
% Save the output of this step

dateStr = datestr(now, 'yyyymmdd_HHMMSS');
save(['Outputs\Step2_Forcasts',num2str(h),'_',dateStr,'.mat'])

%% Calculate MSE from stored forecast error surfaces

MSurf = {MSurf_MP, MSurf_NP, MSurf_FAR1, MSurf_PCA_VAR, MSurf_PCA_KNN, ...
         MSurf_DS_VAR, MSurf_DS_KNN, MSurf_MP_VAR, MSurf_MP_KNN, MSurf_RF};

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

Models = {'MF','NF','FAR', 'PCA VAR','PCA KNN', 'DS VAR','DS KNN', ...
          'MP VAR','MP KNN', 'DS RF'};

mean(MSE_All)

%% Create the figure with the Box plot
% Plot figure in the paper: box plots of forecasts

fig1    = figure(1);
FS      = 14;
canvas  = [100, 100, 1000, 350];
set(fig1, 'Position', canvas, 'PaperPositionMode', 'auto');
bh      = boxplot(MSE_All, 'Labels', Models, 'whisker', 2);
set(bh(7,:), 'Visible', 'off');
ylabel('MSE', 'FontSize', FS);
xtickangle(30)
%Set y-axis limits and font size for labels
ylim([0 250]);

exportgraphics(fig1, ['Outputs/BoxplotsHorizon',num2str(h),'.pdf'], ...
              'BackgroundColor', 'none', 'Resolution', 300);

%% MSE surfaces plots

DTOzone     = delaunayTriangulation([LonOz,LatOz]);
CleanDTOz   = CleanTriangulation(DTOzone,[ConstrReg(:,2),ConstrReg(:,1)],Threshold);
Region      = [ConstrReg(:,2),ConstrReg(:,1)];
RegionBord  = polyshape(Region);
MSESurf     = cellfun(@(M) mean(M.^2,2), MSurf, 'UniformOutput', false);

%--- Compute color axis limits across all models ---
all_mse     = cell2mat(MSESurf');
bottomMSE   = min(all_mse, [], 'all');
topMSE      = max(all_mse, [], 'all');

%--- Plot MSE surfaces ---

figMSE2 = figure(2);

for m = 1:length(Models)

    subplot(2,5,m)
    Season = Grid2Func(MSESurf{m}, DTOzone, CleanDTOz);
    hold on
        plot(RegionBord, 'FaceColor', 'none');
        plot(Season, [], [], [], 100);
    hold off
    axis equal tight;
    colormap(jet);
    clim([bottomMSE topMSE]);
    colorbar;
    view(2);
    title(Models{m});
    xlabel('Easting');
    ylabel('Northing');

end

% 
% %--- Resize and export figure ---
% figLow.Position = [100, 100, 1200, 600];  
% figUp.Position = [100, 100, 1200, 600]; 
% exportgraphics(figLow, ['Outputs/CI_2p5_Surfaces_Horizon', num2str(h), '.pdf'], ...
%     'BackgroundColor', 'none', 'Resolution', 300, 'ContentType', 'vector');
% exportgraphics(figUp, ['Outputs/CI_97p5_Surfaces_Horizon', num2str(h), '.pdf'], ...
%     'BackgroundColor', 'none', 'Resolution', 300, 'ContentType', 'vector');