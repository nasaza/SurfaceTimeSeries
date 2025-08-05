clear all;
clc;

%% Add Libraries, Data and Projections
addpath Data/
addpath AddFunc
load('Data\SeasComp');
load('Data\SeasonAdjData');
ConstrReg = csvread('Data/GeoConstraints/DE_Constraints.csv');




% project data on the plain (Gauss-Krüger Zone 3)
wgs84           = geocrs(4326);
proj            = projcrs(31467);
[LonOz,  LatOz] = projfwd(proj, LatOz, LonOz);
[x, y]          = projfwd(proj, ConstrReg(:,1), ConstrReg(:,2));
ConstrReg       = [y,x];



%% Preparation of the FEM 
Threshold = 0.75;
DTOzone     = delaunayTriangulation([LonOz,LatOz]);
CleanDTOz   = CleanTriangulation(DTOzone,[ConstrReg(:,2),ConstrReg(:,1)],Threshold);


bottom   = min(SeasOz,[],'all');
top      = max(SeasOz,[],'all');
fig      = figure;
months   = {'January', 'February', 'March', 'April', 'May', 'June', ...
          'July', 'August', 'September', 'October', 'November', 'December'};

Region      = [ConstrReg(:,2),ConstrReg(:,1)];
RegionBord  = polyshape(Region);
 

%% Plotting Figure 2

for m = 1:12 
            Season = Grid2Func(SeasOz(:,m),DTOzone,CleanDTOz);
            subplot(4,3,m)
            hold on
            plot(RegionBord,'FaceColor', 'none');
            plot(Season,[],[],[],100);
            % geoshow(LatOz,LonOz,SeasOz(:,m),'DisplayType','surface')
            % demcmap(SeasOz(:,m))
            view(2);
            hold off
            axis equal tight;
            colormap(jet);
            caxis([bottom top])
            colorbar;
            title(months{m});          
            xlabel('Easting');
            ylabel('Northing');
            ax = gca;
            ax.XTick = [];
            ax.YTick = [];
end

fig.Position = [100, 100, 900, 1200]; % [left, bottom, width, height] 

% Adjust the axes to minimize white space
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];


exportgraphics(fig,['Outputs/SeasonalComponets.pdf'],'BackgroundColor','none','Resolution',300,'ContentType', 'vector')   
