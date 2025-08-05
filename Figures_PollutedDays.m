%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Forecasts Comparison
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
clc;

%% Add Libraries
addpath AddFunc
addpath Data

%% Step 1: Read Data

load('Data\SeasonAdjData'); % Data on the grid seasonaly addjusted
load('Data\FTSs');          % Data difen in the Funtional Form (seasonaly addjusted)
ConstrReg = csvread('Data/GeoConstraints/DE_Constraints.csv');
load('Data\SeasComp');      % Loads seasonal componet


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
H           = 165;  %set forecast horizon

[~,T]       = size(OzoneS);
K_max       = 15;
p_max       = 5;    % Maximum nuber of lag considered in the analysis
KNN_max     = 100;   % Maximum number of KNN in the calibraiton

FEMBasOz    = OzoneBasis;


L_set = [3 1 1 1 1]; %Ozone, Sun, Precip, Temp, Wind
m_set = [2 1 1 1 1];
 
   
DTOzone     = delaunayTriangulation([LonOz,LatOz]);
CleanDTOz   = CleanTriangulation(DTOzone,[ConstrReg(:,2),ConstrReg(:,1)],0.75);
Region      = [ConstrReg(:,2),ConstrReg(:,1)];
RegionBord  = polyshape(Region);


        
%% Step 3: Forecasting with Different Methods
h=1;
tic
for i=1:H 
    if (i==36)|(i==38)

        %runs loop with updating information 
        % information can be updates with all parameters or without updating
        % K,p,L,q
    %% Recover part of the data to be used for estimation
        nobs_i      = T-H+i-1;
        TrueVal     = OzoneS(:,T-H+i)+SeasOz(:,8);

        OzoneFTSi   = fd(OzoneCoef(:,1:T-H+i-1),OzoneBasis);

    %% Functional Persepctive (FP): fPCA for the available observations T-(H-1)+i

        pcastrOz_i  = pca3D({OzoneCoef(:,1:T-H+i-1),OzoneBasis}, K_max, 1); 
        pcastrS_i   = pca3D({SunCoef(:,1:T-H+i-1),SunBasis}, K_max, 1); 
        pcastrP_i   = pca3D({PrecipCoef(:,1:T-H+i-1),PrecipBasis}, K_max, 1); 
        pcastrT_i   = pca3D({TempCoef(:,1:T-H+i-1),TempBasis}, K_max, 1); 
        pcastrW_i   = pca3D({WindCoef(:,1:T-H+i-1),WindBasis}, K_max, 1); 
        PCA_addreg  = {pcastrS_i,pcastrP_i,pcastrT_i,pcastrW_i};

    %% Method 1: (FP) Mean Predictor

        MSurf1      = mean(OzoneS(:,1:T-H+i-1),2)+SeasOz(:,8);

    %% Method 2: (FP) Naive Predictor

        MSurf2      = OzoneS(:,T-H+i-1)+SeasOz(:,8);

    %% Method 3: (FP) FAR(1)

        FuncPredCoeff   = SFAR1(pcastrOz_i,5,h);
        FuncPred        = fd(FuncPredCoeff,FEMBasOz);
        FPointEval      = eval_FEM_fd(LonOz,LatOz,FuncPred);
        MSurf3          = FPointEval+SeasOz(:,8);

    %% Method 4: (FP) Linear Forecast with scree plot estimates

        FuncPredCoeff   = FFM_VARX(pcastrOz_i,PCA_addreg,L_set,m_set,h);
        FuncPred        = fd(FuncPredCoeff,FEMBasOz);
        FPointEval      = eval_FEM_fd(LonOz,LatOz,FuncPred);
        MSurf4          = FPointEval+SeasOz(:,8);

    %% Method 5: (Multivariate) Linear Forecast with scree plot estimates

        FPredGrid     = DFFM_VARX_GRID(OzoneS(:,1:T-H+i-1),WindS(:,1:T-H+i-1),...
                            SunS(:,1:T-H+i-1),PrecipS(:,1:T-H+i-1),...
                            TempS(:,1:T-H+i-1),L_set,m_set,h);
       MSurf5         = FPredGrid+SeasOz(:,8);


    % %% Method 6: (FP) Nonlinear Forecast with KNN
    % % The code for model 6 has two parts. One with searching over optimal KNNs
    % % post selection and one runs with KNN=10
    % %-> Part 1: Estimating number of neigbours
    %     K_NNv        = zeros(KNN_max,1) ;
    %     pcastrKNN_i  = pca3D({OzoneCoef(:,1:T-H+i-2),OzoneBasis}, K_max, 1); 
    %     for knn = 1:KNN_max
    %         FuncPredCoeff  = FFM_KNN(pcastrKNN_i,L_set(1),m_set(1),knn);
    %         FuncPred       = fd(FuncPredCoeff,FEMBasOz);
    %         FPointEval     = eval_FEM_fd(LonOz,LatOz,FuncPred);
    %         K_NNv(knn,1)   = mean((OzoneS(:,T-H+i-1)-FPointEval).^2);
    %     end
    %     [~,K_min]     = min(K_NNv);
    %     KNNs(i,:)     = K_min;
    % %-> Part 2    
    %     FuncPredCoeff  = FFM_KNN(pcastrOz_i,L_set(1),m_set(1),K_min);
    %     FuncPred       = fd(FuncPredCoeff,FEMBasOz);
    %     FPointEval     = eval_FEM_fd(LonOz,LatOz,FuncPred);
    %     MSurf6         = FPointEval+SeasOz(:,8);


    %% Method 7: (MP) Nonlinear Forecast with KNN 
    % The code for model 7 (as for model 6) has two parts. One with searching over optimal KNNs
    % post selection and one runs with KNN=10
    %-> Part 1
    %     K_NNv       = zeros(KNN_max,1);
    %     pcastr_i    = FactorDecompFTSGrid(OzoneS(:,1:T-H+i-2), K_max, 1);
    %     for knn = 1:KNN_max
    %         FPredGrid      = DFFM_KNN_GRID(pcastr_i,L_set(1),m_set(1),knn);
    %         K_NNv(knn,1)   = mean((OzoneS(:,T-H+i-1)-FPredGrid).^2);
    %     end
    %     [~,K_min]      = min(K_NNv);
    %     KNNGs(i,:)     = K_min;
    % %-> Part 2
    %     pcastr_i       = FactorDecompFTSGrid(OzoneS(:,1:T-H+i-1), K_max, 1);
    %     FPredGrid      = DFFM_KNN_GRID(pcastr_i,L_set(1),m_set(1),K_min);
    %     MSurf7         = FPredGrid+SeasOz(:,8);   
        
     %% Plot forecast surfaces     
        fig = figure;
%         Models      = {'True Surf','MP', 'NP', 'FAR', 'LGF', 'LSF', 'NGF', 'NSF'};
%         SurfArray   = {TrueVal, MSurf1, MSurf2, MSurf3, MSurf5, MSurf4, MSurf7, MSurf6};
%         
%         bottom   = min([TrueVal, MSurf1, MSurf2, MSurf3, MSurf5, MSurf4, MSurf7, MSurf6],[],'all');
%         top      = max([TrueVal, MSurf1, MSurf2, MSurf3, MSurf5, MSurf4, MSurf7, MSurf6],[],'all');

        Models      = {'True Surf', 'FAR', 'LGF', 'LSF'};
        SurfArray   = {TrueVal, MSurf3, MSurf5, MSurf4};
        
        bottom   = min([TrueVal, MSurf3, MSurf5, MSurf4],[],'all');
        top      = max([TrueVal, MSurf3, MSurf5, MSurf4],[],'all');

        for m = 1:length(Models) 
                    Season  = Grid2Func(SurfArray{m},DTOzone,CleanDTOz);
                    subplot(1,4,m)
%                     subplot(1,4,m)
                    hold on
                    plot(RegionBord,'FaceColor', 'none');
                    plot(Season,[],[],[],100);
                    hold off
                    axis equal tight
                    colorbar
                    colormap(jet)    
                    % if m==length(Models) 
                    caxis([min(SurfArray{m},[],'all') max(SurfArray{m},[],'all')])
                    colorbar;
                    % end
                    view(2);
                    title(Models{m});
                    xlabel('Easting');
                    ylabel('Northing');
                    % ax = gca;
                    % ax.XTick = [];
                    % ax.YTick = [];
        end
        fig.Position = [100, 100, 1200, 300]; % [left, bottom, width, height] 
        exportgraphics(fig,['Outputs/FSurfaces',num2str(i),'.pdf'],'BackgroundColor','none','Resolution',300,'ContentType', 'vector')   

    end
    i
end
toc


%% MSE surfaces
