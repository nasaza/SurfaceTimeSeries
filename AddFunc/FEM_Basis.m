 function FEMbasis = FEM_Basis(DT,CleanedDT)
    % Create FEM basis based on DT trinagulation with cleaned DT 
    
    Coord    = [DT.Points(:,1), DT.Points(:,2)];    % Points/Coordiates
    TriangD  = DT.ConnectivityList;                 % Triangles
    TriangD  = TriangD(logical(CleanedDT),:);       % Cleaned Triangles
    FEMbasis = create_FEM_basis(Coord, [], TriangD, 2);
 
end