clear all;
clc;

%% Short Info:
% This file creates functiona time sereis

%% Add Libraries and Data
scriptFolder = fileparts(mfilename('fullpath')); % Repository path
rootFolder   = scriptFolder;

addpath(fullfile(rootFolder,'AddFunc'));
dataFolder   = fullfile(rootFolder,'Data');
outFolder    = fullfile(rootFolder,'Outputs');
if ~exist(outFolder,'dir')
    mkdir(outFolder);
end

%% Step 1: Read Data and Project
% This step loads seasonally adjusted ozone concentration data observed
% at different stations, along with other weather variables used in the analysis
% and geographic borders of Germany

load(fullfile(dataFolder,'SeasonAdjData.mat'));  % seasonally adjusted data on observed locations
ConstrReg = readmatrix(fullfile(dataFolder,'GeoConstraints','DE_Constraints.csv'));


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

%% Step 2: Triangulation of the selected geographic area
% This step defines the geographical area for our analysis. In particular, 
% it creates a triangulation of the area for each variable and generates Figure A.1 
% in the paper's appendix.

DTOzone     = delaunayTriangulation([LonOz,LatOz]);
CleanDTOz   = CleanTriangulation(DTOzone,[ConstrReg(:,2),ConstrReg(:,1)],Threshold);
% [RemP, regP] = TriangStats(DTOzone, [ConstrReg(:,2),ConstrReg(:,1)], Threshold);

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

fig1 = figure(1);
plot(RegionBord,'FaceColor', 'none'); 
hold on;
triplot(DTOzone(logical(CleanDTOz),:),DTOzone.Points(:,1), DTOzone.Points(:,2));
axis equal tight;
hold off;    
xlabel('Easting (meters)');
ylabel('Northing (meters)');
exportgraphics(fig1, fullfile(outFolder,'FigureAppTriangulation.pdf'), ...
    'BackgroundColor','none', 'Resolution',300);

%% Step 3: Create Functinal Data  
% Creates surface/functional observations from gridded data. It uses
% the "Grid2Func" function from the "AddFunc" folder, which is included
% in this package.
% This step may take some time. To speed up the procedure, once surface
% observations are created for the first time, the required coefficient
% matrices and basis objects are saved in Data/FTSs.mat.
% To reload the generated objects later, use: load(fullfile(dataFolder,'FTSs.mat'));

[~,OzoneCoef,OzoneBasis]      = Grid2Func(OzoneS,DTOzone,CleanDTOz);
[~,SunCoef,SunBasis]          = Grid2Func(SunS,DTSun,CleanDTSun);
[~,PrecipCoef,PrecipBasis]    = Grid2Func(PrecipS,DTPrecip,CleanDTPrec);
[~,TempCoef,TempBasis]        = Grid2Func(TempS,DTTemp,CleanDTTemp);
[~,WindCoef,WindBasis]        = Grid2Func(WindS,DTWind,CleanDTWind);
save(fullfile(dataFolder,'FTSs.mat'),'OzoneCoef','OzoneBasis','SunCoef',...
    'SunBasis','PrecipCoef','PrecipBasis','TempCoef','TempBasis',...
    'WindCoef','WindBasis');
