clear all;
clc;

%% Short Info:
% This file creates functiona time sereis

%% Add Libraries and Data
addpath AddFunc
addpath Data

%% Step 1: Read Data
% This step loads seasonally adjusted ozone concentration data observed
% at different stations, along with other weather variables used in the analysis
% and geographic borders of Germany

load('Data\SeasonAdjData');
ConstrReg = csvread('Data/GeoConstraints/DE_Constraints.csv');
mkdir('Outputs');  % creates a folder to save all of the outputs 

%% Step 2: Triangulation of the selected geographic area
% This step defines the geographical area for our analysis. In particular, 
% it creates a triangulation of the area for each variable and generates Figure 7 in the paper's appendix.

DTOzone     = delaunayTriangulation([LonOz,LatOz]);
OutDTs      = [2,5,20,35,69,70,76,80,81,128,131,133,134,135,159,...
                        173,239,248,261,308,312,315,322,323,324];
CleanDTOz           = ones(326,1);
CleanDTOz(OutDTs)   = zeros(length(OutDTs),1);   


DTSun       = delaunayTriangulation([XSun,YSun]);
CleanDTSun  = CleanTriangulation(DTSun,[ConstrReg(:,2),ConstrReg(:,1)]);

DTPrecip    = delaunayTriangulation([XPrecip,YPrecip]);
CleanDTPrec = CleanTriangulation(DTPrecip,[ConstrReg(:,2),ConstrReg(:,1)]);

DTTemp      = delaunayTriangulation([XTemp,YTemp]);
CleanDTTemp = CleanTriangulation(DTTemp,[ConstrReg(:,2),ConstrReg(:,1)]);

DTWind      = delaunayTriangulation([XWind,YWind]);
CleanDTWind = CleanTriangulation(DTWind,[ConstrReg(:,2),ConstrReg(:,1)]);

Region      = [ConstrReg(:,2),ConstrReg(:,1)];
RegionBord  = polyshape(Region);

figure(1)
plot(RegionBord)    
hold on;
triplot(DTOzone(logical(CleanDTOz),:),DTOzone.Points(:,1), DTOzone.Points(:,2));
axis equal;
hold off;    

%% Step 3: Create Functinal Data  
% Creates surface/functional observations from gridded data. It uses
% the "Grid2Func" function from the "AddFunc" folder, which is included
% in this package.
% This step may take some time. To speed up the procedure, once surface
% observations are created for the first time, they are saved in the
% "Output" folder and can be reloaded on demand. To skip Step 3 and
% load previously generated data instead, comment out this block and
% uncomment the last line:  load('Data\FTSs');.

[OzoneFTS,OzoneCoef,OzoneBasis]    = Grid2Func(OzoneS,DTOzone,CleanDTOz);
[SunFTS,SunCoef,SunBasis]          = Grid2Func(SunS,DTSun,CleanDTSun);
[PrecipFTS,PrecipCoef,PrecipBasis] = Grid2Func(PrecipS,DTPrecip,CleanDTPrec);
[TempFTS,TempCoef,TempBasis]       = Grid2Func(TempS,DTTemp,CleanDTTemp);
[WindFTS,WindCoef,WindBasis]       = Grid2Func(WindS,DTWind,CleanDTWind);
OzoneFTSobj                        = {OzoneFTS,OzoneCoef,OzoneBasis};
SunFTSobj                          = {SunFTS,SunCoef,SunBasis};
PrecipFTSobj                       = {PrecipFTS,PrecipCoef,PrecipBasis};
TempFTSobj                         = {TempFTS,TempCoef,TempBasis};
WindFTSobj                         = {WindFTS,WindCoef,WindBasis};
save('Data\FTSs','OzoneFTS','OzoneCoef','OzoneBasis','OzoneFTSobj',...
    'SunFTS','SunCoef','SunBasis','SunFTSobj','PrecipFTS',...
    'PrecipCoef','PrecipBasis','PrecipFTSobj','TempFTS',...
    'TempCoef','TempBasis','TempFTSobj','WindFTS',...
    'WindCoef','WindBasis','WindFTSobj');

% If FTS are created than we can load them quicker as 
 % load('Data\FTSs');          % Data difen in the Funtional Form (seasonaly addjusted)


%% Additional Plots 
% This block obtains additional information about the characteristics of our surface data.
% Specifically, it plots scree plots for each variable to help determine how many PCA components to include in our analysis.
% Replicates Figure 9 in the paper's appendix.

K_max       = 15;    
pcastrOzone = pca3D(OzoneFTS, K_max, 1);     
pcastrTemp  = pca3D(TempFTS(1:200), K_max, 1);     
pcastrSun   = pca3D(SunFTS(1:200), K_max, 1);     
pcastrWind  = pca3D(WindFTS(1:200), K_max, 1);     
pcastrPrecip= pca3D(PrecipFTS(1:200), K_max, 1);     

% plot a scree plot of explaind variance for an intuition
h = findobj('type','figure');
n = length(h);

% scree plots

figure(1);
subplot(3,2,1)
    plot(pcastrOzone.varprop(1:K_max));
    title('Ozone');
subplot(3,2,2)
    plot(pcastrTemp.varprop(1:K_max));
    title('Temp');
subplot(3,2,3)
    plot(pcastrSun.varprop(1:K_max));
    title('Sun');
subplot(3,2,4)
    plot(pcastrWind.varprop(1:K_max));
    title('Wind');
subplot(3,2,5)
    plot(pcastrPrecip.varprop(1:K_max))
    title('Precip');


    % %% Lag and dimension selection based on Aue 
    % 
    % [~,ANH] = InformCriteria(pcastrOzone,15,10);


fh=figure(2);
subplot(2,2,1)
    autocorr(pcastrOzone.pcascr(:,1));
    title('1st Score Series');
subplot(2,2,2)
    autocorr(pcastrOzone.pcascr(:,2));
    title('2nd Score Series');
subplot(2,2,3)
    autocorr(pcastrOzone.pcascr(:,3));
    title('3rd Score Series');
subplot(2,2,4)
    autocorr(pcastrOzone.pcascr(:,4));
    title('4th Score Series');

exportgraphics(fh, ['Outputs/AutocorrScores.pdf'], 'BackgroundColor', 'none', 'Resolution', 300);
