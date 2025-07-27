
%% Add Libraries and Data
addpath Data/
addpath AddFunc
load('Data\SeasComp');
ConstrReg = csvread('Data/GeoConstraints/DE_Constraints.csv');

%% Preparation of the FEM and data
DTOzone     = delaunayTriangulation([LonOz,LatOz]);
CleanDTOz   = ones(326,1);
OutDTs      = [2,5,20,35,69,70,76,80,81,128,131,133,134,135,159,...
                        173,239,248,261,308,312,315,322,323,324];

CleanDTOz(OutDTs)   = zeros(length(OutDTs),1);     


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
            colormap(jet);
            caxis([bottom top])
            colorbar;
            title(months{m});
            xlabel('');
            ylabel('');
            ax = gca;
            ax.XTick = [];
            ax.YTick = [];
end

fig.Position = [100, 100, 600, 500]; % [left, bottom, width, height] 

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
