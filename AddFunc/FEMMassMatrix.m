function K0mat = FEMMassMatrix(basisobj)
% FEMMassMatrix
% Extracts/recomputes the FEM mass matrix for an FDA FEM basis object.
%
% K0mat(i,j) = int_D phi_i(s) phi_j(s) ds
%
% This reuses the FDA package's internal mass(nodeStruct) routine.

if ~strcmp(getbasistype(basisobj), 'FEM')
    error('Argument basisobj is not of type FEM.');
end

params = getbasispar(basisobj);

nodeStruct.order     = params.order;
nodeStruct.nodes     = params.nodes;
nodeStruct.nodeindex = params.nodeindex;
nodeStruct.J         = params.J;
nodeStruct.metric    = params.metric;

K0mat = mass(nodeStruct);

% Numerical symmetrization
K0mat = (K0mat + K0mat')/2;

end