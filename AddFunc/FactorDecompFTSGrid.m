function pcastr = FactorDecompFTSGrid(FTS_grid, npca, centerfd)
%  Decomposes surface FTS evalauted on the sampled geo-locations into factors 
% Inputs:
%  - FTS_grid       functional time series evaluated on the geo "grid" given  
%                     as JxT matrix, with J number of sampling points and T sample size 
%  - npca           number of pca components/factors to be computed;
%  - centerfd       centered or not centered data; 0 for not centered and 1
%                   for centered


%% Step 1: De meaning of the FTS
[~,T]    = size(FTS_grid);
    
if nargin < 3
    centerfn = 1;   %  subtract mean from data before PCA
end

if centerfd == 1 
    meanfd   = mean(FTS_grid,2);
    FTS_grid = FTS_grid-kron(ones(1,T), meanfd);
end

%% Step 2: Calculate estimate of the Cov of FTS_grid

CovX       = 1/T*(FTS_grid'*FTS_grid);

%% Step 3: calculate  nharm eigenvectors of CovX

[Theta,D]  = eig(CovX); % Theta is a matrix eigenvectors
                        % D if diagonal matrix with eigenvalues
eigvals    = diag(D);
[~,inds]   = sort(eigvals);
eigvals    = eigvals(flip(inds));
Theta      = Theta(:,flip(inds));
Theta      = Theta(:,1:npca);

%% Saving outputs 
scores     = sqrt(T)*Theta;
varprop    = 1/sum(eigvals)*cumsum(eigvals);
EigenVec   = (1/T)*FTS_grid*scores;

pcastr.values  = eigvals(1:npca);
pcastr.scores  = scores(:,1:npca);
pcastr.varprop = varprop(1:npca);
pcastr.EigVec  = EigenVec(:,1:npca);
pcastr.mean    = meanfd;
