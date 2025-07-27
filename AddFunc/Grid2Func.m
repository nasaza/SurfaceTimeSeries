function [FTSobj,coefs,FEMbasis]=Grid2Func(EnvData,DT,CleanedDT)
    % Seeting up initial values  
    
    Coord   = [DT.Points(:,1), DT.Points(:,2)];    % Points/Coordiates
    TriangD = DT.ConnectivityList;                 % Triangles
    TriangD = TriangD(logical(CleanedDT),:);       % Cleaned Triangles
    np      = size(Coord,1);                       % Number of Points/Coordinates
    nt      = size(TriangD,1);                     % Number of triangles
    T       = size(EnvData,2);                     % Number of observations

    % Set up the FEM basis object 
    order       = 2;                  % Order of local polynomials
    FEMbasisobj = create_FEM_basis(Coord, [], TriangD, order);
    Envir_fd    = fd(zeros(getnbasis(FEMbasisobj),1),FEMbasisobj);  % Setting up an empty FEM functional data object
    FEM_cell    = cell(T,1);                                        % Storage for functional data objects of each daily surface (i.e. spatial smooth over the whole domain for each day) in a cell array

    % -------------------------------------------------------------------------
    % Run spatial smoothing and store smoothed fd object ozone_fd_cell

    lambda   = 0.01;
    coef_FEM = zeros(getnbasis(FEMbasisobj),T);
    for iday=1:T
        data_i              = [(1:np)', EnvData(:,iday)];
        [ozone_fd_i]        = smooth_FEM_fd(data_i,Envir_fd,lambda);        
        FEM_cell{iday,1}    = ozone_fd_i;
        coef_FEM(:,iday)    = getcoef(FEM_cell{iday,1});
    end
    FTS = fd(coef_FEM,FEMbasisobj);
    
    if nargout == 1
        FTSobj  = FTS;
    elseif nargout == 2
        FTSobj  = FTS;
        coefs   = coef_FEM;
    elseif nargout == 3        
        FTSobj  = FTS;
        coefs   = coef_FEM;
        FEMbasis= FEMbasisobj;
    end
end