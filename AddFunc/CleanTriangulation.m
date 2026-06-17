function TrigInd = CleanTriangulation(DT,Constarains,threshold)
% This function cleans Delaunay Triangulation, DT, from triangles outside of
% geographical boarders pecified in Contrains
% 
% Inputs:
%       DT - Delaunay Triangulation
%       Constrains - (n x 2) matrix representing contraining polygon; first
%       column is Lon-coordinate and the second one is Lat-coordinate

if nargin < 3
    threshold = 1; % default value
end


% create Region boarders polyshape which contarins triangulation
RegionB = polyshape(Constarains(:,1),Constarains(:,2)); 

% Extrat each of the triangles and check if it is inside of region boarder RegionB

[TrigN,~] = size(DT);
TrigInd   = zeros(TrigN,1);
for ii = 1:TrigN
    Trig_i      = DT.Points(DT(ii,:),:);
    Trig_i      = [Trig_i;Trig_i(1,:)]; % Add the first coordinate to create closed polygon
    PolyTrig    = polyshape(Trig_i);

    % Compute intersection area
    PolyIntst   = intersect(RegionB, PolyTrig);
    AreaTrig    = area(PolyTrig);
    AreaIntersect = area(PolyIntst);
    
    % Keep triangle if most of it is inside the region
    TrigInd(ii) = AreaIntersect / AreaTrig >= threshold;
end

end