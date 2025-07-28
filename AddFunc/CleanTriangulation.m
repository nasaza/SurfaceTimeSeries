function TrigInd = CleanTriangulation(DT,Constarains)
% This function cleans Delaunay Triangulation, DT, from triangles outside of
% geographical boarders pecified in Contrains
% 
% Inputs:
%       DT - Delaunay Triangulation
%       Constrains - (n x 2) matrix representing contraining polygon; first
%       column is Lon-coordinate and the second one is Lat-coordinate

%See TestingCode2 for an Examples

% create Region boarders polyshape which contarins triangulation
RegionB = polyshape(Constarains(:,1),Constarains(:,2)); 

% Extrat each of the triangles and check if it is inside of region boarder RegionB

[TrigN,~]=size(DT);
TrigInd  = zeros(TrigN,1);
for ii = 1:TrigN
    Trig_i      = DT.Points(DT(ii,:),:);
    Trig_i      = [Trig_i;Trig_i(1,:)]; % Add the first coordinate to create closed polygon
    PolyTrig    = polyshape(Trig_i);
    Z1          = subtract(RegionB,PolyTrig);
    Z2          = subtract(PolyTrig,RegionB);        
    TrigInd(ii)= (Z1.NumRegions~=0)&(Z2.NumRegions==0);
end
% 
% CleanedDT = DT;
% CleanedDT.ConnectivityList = [];

end