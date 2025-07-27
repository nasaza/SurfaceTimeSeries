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

FEMBasOz    = OzoneBasis;

MSE_MP      = zeros(N_Fors,1);
MSE_NP      = zeros(N_Fors,1);
MSE_ANH     = zeros(N_Fors,1);
MSE_LF_FP   = zeros(N_Fors,1);
MSE_OS_HQ   = zeros(N_Fors,1);
MSE_ANHR    = zeros(N_Fors,1);
MSE_OSR     = zeros(N_Fors,1);
MSE_OS_HQR  = zeros(N_Fors,1);
MSE_HJ      = zeros(N_Fors,1);
MSE_VARX_g  = zeros(N_Fors,1);
MSE_KNN     = zeros(N_Fors,1);
MSE_KNNG    = zeros(N_Fors,1);
MSE_FAR1    = zeros(N_Fors,1);
MSE_RF      = zeros(N_Fors,1);


MSurf1      = zeros(size(OzoneS,1),1); 
MSurf2      = zeros(size(OzoneS,1),1); 
MSurf3      = zeros(size(OzoneS,1),1); 
MSurf4      = zeros(size(OzoneS,1),1); 
MSurf5      = zeros(size(OzoneS,1),1); 
MSurf6      = zeros(size(OzoneS,1),1); 
MSurf7      = zeros(size(OzoneS,1),1); 
MSurf8      = zeros(size(OzoneS,1),1); 

% Run Aue et al (2015) to select the follwoing values (function InformCriteria(PCA_Str,K_max,p_max))
L_set = [3 1 1 1 1]; % Ozone, Sun, Precip, Temp, Wind
m_set = [2 1 1 1 1];


%% Step 3: Forecasting with Different Methods
% Calculates forecasts for horizon h over N_Fors iterations, updating the sample each time
% and saving the forecast outputs in the corresponding variables

tic

for i=1:(H+1-h)    
%% Part of the data to be used for estimation
    nobs_i      = T-H+i-1;
    TrueVal     = OzoneS(:,nobs_i+h);    
    OzoneFTSi   = fd(OzoneCoef(:,1:nobs_i),OzoneBasis);
    
%% Functional Persepctive (FP): fPCA for the available observations T-(H-1)+i

    pcastrOz_i  = pca3D({OzoneCoef(:,1:nobs_i),OzoneBasis}, K_max, 1); 
    pcastrS_i   = pca3D({SunCoef(:,1:nobs_i),SunBasis}, K_max, 1); 
    pcastrP_i   = pca3D({PrecipCoef(:,1:nobs_i),PrecipBasis}, K_max, 1); 
    pcastrT_i   = pca3D({TempCoef(:,1:nobs_i),TempBasis}, K_max, 1); 
    pcastrW_i   = pca3D({WindCoef(:,1:nobs_i),WindBasis}, K_max, 1); 
    PCA_addreg  = {pcastrS_i,pcastrP_i,pcastrT_i,pcastrW_i};
    
%% Method 1: (FP) Mean Predictor

    MSE_MP(i)   = mean((TrueVal-mean(OzoneS(:,1:nobs_i),2)).^2);
    MSurf1      = MSurf1+((TrueVal-mean(OzoneS(:,1:nobs_i),2)).^2)/N_Fors;

%% Method 2: (FP) Naive Predictor

    MSE_NP(i)   = mean((TrueVal-OzoneS(:,nobs_i)).^2);
    MSurf2      = MSurf2+((TrueVal-OzoneS(:,nobs_i)).^2)/N_Fors;
    
    
%% Method 3: (FP) FAR(1)

    FuncPredCoeff   = SFAR1(pcastrOz_i,5,h);
    FuncPred        = fd(FuncPredCoeff,FEMBasOz);
    FPointEval      = eval_FEM_fd(LonOz,LatOz,FuncPred);
    MSE_FAR1(i)     = mean((TrueVal-FPointEval).^2); 
    MSurf3          = MSurf3+((TrueVal-FPointEval).^2)/N_Fors;
    
%% Method 4: (FP) Linear Forecast with scree plot estimates

    FuncPredCoeff   = FFM_VARX(pcastrOz_i,PCA_addreg,L_set,m_set,h);
    FuncPred        = fd(FuncPredCoeff,FEMBasOz);
    FPointEval      = eval_FEM_fd(LonOz,LatOz,FuncPred);
    MSE_LF_FP(i)    = mean((TrueVal-FPointEval).^2);
    MSurf4          = MSurf4+((TrueVal-FPointEval).^2)/N_Fors;

%% Method 5: (Multivariate) Linear Forecast with scree plot estimates

    FPredGrid     = DFFM_VARX_GRID(OzoneS(:,1:T-H+i-1),WindS(:,1:T-H+i-1),...
                        SunS(:,1:T-H+i-1),PrecipS(:,1:T-H+i-1),...
                        TempS(:,1:T-H+i-1),L_set,m_set,h);
    MSE_VARX_g(i) = mean((TrueVal-FPredGrid).^2);
    MSurf5        = MSurf5+((TrueVal-FPredGrid).^2)/N_Fors;

     
%% Method 6: (FP) Nonlinear Forecast with KNN
%-> Part 1: Estimating number of neigbours

K_NNv        = zeros(KNN_max,1) ;
    pcastrKNN_i  = pca3D({OzoneCoef(:,1:nobs_i-h),OzoneBasis}, K_max, 1); 
    for knn = 1:KNN_max
        FuncPredCoeff  = FFM_KNN(pcastrKNN_i,L_set(1),m_set(1),knn,h);
        FuncPred       = fd(FuncPredCoeff,FEMBasOz);
        FPointEval     = eval_FEM_fd(LonOz,LatOz,FuncPred);
        K_NNv(knn,1)   = mean((OzoneS(:,nobs_i)-FPointEval).^2);
    end
    [~,K_min]     = min(K_NNv);
    KNNs(i,:)     = K_min;
%-> Part 2    
    FuncPredCoeff  = FFM_KNN(pcastrOz_i,L_set(1),m_set(1),K_min,h);
    FuncPred       = fd(FuncPredCoeff,FEMBasOz);
    FPointEval     = eval_FEM_fd(LonOz,LatOz,FuncPred);
    MSE_KNN(i)     = mean((TrueVal-FPointEval).^2);
    MSurf6         = MSurf6+((TrueVal-FPointEval).^2)/N_Fors;
     
    
%% Method 7: (MP) Nonlinear Forecast with KNN 
% to be reconstructed
%-> Part 1
    K_NNv       = zeros(KNN_max,1);
    pcastr_i    = FactorDecompFTSGrid(OzoneS(:,1:T-H+i-2), K_max, 1);
    for knn = 1:KNN_max
        FPredGrid      = DFFM_KNN_GRID(pcastr_i,L_set(1),m_set(1),knn);
        K_NNv(knn,1)   = mean((OzoneS(:,T-H+i-1)-FPredGrid).^2);
    end
    [~,K_min]      = min(K_NNv);
    KNNGs(i,:)     = K_min;
% %-> Part 2
    pcastr_i       = FactorDecompFTSGrid(OzoneS(:,1:T-H+i-1), K_max, 1);
    FPredGrid      = DFFM_KNN_GRID(pcastr_i,L_set(1),m_set(1),K_min);
    MSE_KNNG(i)    = mean((TrueVal-FPredGrid).^2);
    MSurf7         = MSurf7+((TrueVal-FPredGrid).^2)/N_Fors;      

%% Method 8: (Multivariate) Linear Forecast with Random Forest

    FuncPredCoeff   = FRandomForest(pcastrOz_i,L_set(1),h);
    FuncPred        = fd(FuncPredCoeff,FEMBasOz);
    FPointEval      = eval_FEM_fd(LonOz,LatOz,FuncPred);
    MSE_RF(i)       = mean((TrueVal-FPointEval).^2);
    MSurf8          = MSurf8+((TrueVal-FPointEval).^2)/N_Fors;
    i
end
toc

%% Analyzing the forecasting performance
% This step compares the performance of different methods based on their MSE 
% and generates Figures 4 and 5 in the main text of the paper

% Get the current date and time as a string
dateStr = datestr(now, 'yyyymmdd_HHMMSS');
save(['Outputs\Step2_Forcasts',num2str(h),'_',dateStr,'.mat'])
 
mean([MSE_MP,MSE_NP,MSE_FAR1,MSE_LF_FP,MSE_VARX_g,MSE_KNN,MSE_KNNG,MSE_RF])



% Define font sizes and canvas size
FS = 14;
canvas = [100, 100, 850, 350]; % Adjusted canvas size

% Create the figure
fh = figure(1);
set(fh, 'Position', canvas, 'PaperPositionMode', 'auto');

% Plot the boxplot
bh = boxplot([MSE_MP, MSE_NP, MSE_FAR1, MSE_VARX_g, MSE_LF_FP, MSE_KNNG, MSE_KNN,MSE_RF], ...
            'Labels', {'MP', 'NP', 'FAR', 'LGF', 'LSF', 'NGF', 'NSF','RandFor'}, 'whisker', 2);

% Set y-axis limits and font size for labels
ylim([0 450]);
set(bh(7,:), 'Visible', 'off');
ylabel('MSE', 'FontSize', FS);

exportgraphics(fh, ['Outputs/BoxplotsHorizon',num2str(h),'.pdf'], 'BackgroundColor', 'none', 'Resolution', 300);

%% MSE surfaces
fig2 = figure(2);
DTOzone     = delaunayTriangulation([LonOz,LatOz]);
CleanDTOz   = ones(326,1);
OutDTs      = [2,5,20,35,69,70,76,80,81,128,131,133,134,135,159,...
                        173,239,248,261,308,312,315,322,323,324];
CleanDTOz(OutDTs)   = zeros(length(OutDTs),1);     

% Msurf    = {MSurf1,MSurf2,MSurf3,MSurf4,MSurf5,MSurf6,MSurf7};
Msurf    = {MSurf2,MSurf3,MSurf5,MSurf4,MSurf7,MSurf6};

bottom   = min([MSurf1,MSurf2,MSurf3,MSurf5,MSurf4,MSurf7,MSurf6],[],'all')
top      = max([MSurf1,MSurf2,MSurf3,MSurf4,MSurf5,MSurf6,MSurf7],[],'all')
% Models   = {'MP', 'NP', 'FAR', 'M1', 'M2', 'M3', 'M4'};
Models   = {'NP', 'FAR', 'LGF', 'LSF', 'NGF', 'NSF'};
Region      = [ConstrReg(:,2),ConstrReg(:,1)];
RegionBord  = polyshape(Region);


for m = 1:6 
            Season  = Grid2Func(Msurf{m},DTOzone,CleanDTOz);
            subplot(2,3,m)
            hold on
            plot(RegionBord,'FaceColor', 'none');
            plot(Season,[],[],[],100);
            hold off
            colormap(jet);
            caxis([bottom top])
            colorbar;
            view(2);
            title(Models{m});
            xlabel('');
            ylabel('');
            ax = gca;
            ax.XTick = [];
            ax.YTick = [];
end

fig2.Position = [100, 100, 600, 300]; % [left, bottom, width, height] 
exportgraphics(fig2,['Outputs/MSEsurfacesHorizon',num2str(h),'.pdf'],'BackgroundColor','none','Resolution',300,'ContentType', 'vector')   