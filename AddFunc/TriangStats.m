function [removalRate, regionCoverage] = TriangStats(DT, Constarains, minCoverage)
    RegionPoly = polyshape(Constarains(:,1), Constarains(:,2));
    RegionArea = area(RegionPoly);

    numTriangles = size(DT.ConnectivityList, 1);
    TrigInd      = false(numTriangles, 1);  % preallocate
    coveredArea  = 0;

    for ii = 1:numTriangles
        nodes = DT.ConnectivityList(ii,:);
        coords = DT.Points(nodes,:);
        coords = [coords; coords(1,:)];  % close the triangle
        triPoly = polyshape(coords);

        intersection = intersect(triPoly, RegionPoly);
        A_trig = area(triPoly);
        A_int  = area(intersection);

        if A_trig > 0 && A_int / A_trig >= minCoverage
            TrigInd(ii) = true;
            coveredArea = coveredArea + A_int;
        end
    end

    removalRate    = 100 * (1 - sum(TrigInd)/numTriangles);
    regionCoverage = 100 * (coveredArea / RegionArea);
end
