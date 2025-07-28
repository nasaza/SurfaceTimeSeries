function FuncPredCoeff = SFAR1(PCA_Str,L,h_horizon)
    %% remove the mean before handaling the data!
    % Note PCA_Str has to be recovered using pca3D
    % initial values
    Fscores         = PCA_Str.pcascr;
    Feigval         = PCA_Str.values;
    [nobs,~]        = size(Fscores);
        
    % Calculating kernel matrix 
    A               = Fscores(1:end-1,1:L);
    B               = Fscores(2:end,1:L);
    C               = (((Feigval(1:L)).^(-1)).* (Fscores(end,1:L))')/(nobs-1);

    if ~isscalar(h_horizon) || h_horizon < 1 || h_horizon ~= floor(h_horizon)
        error('Input must be an integer starting from 1.');
    end

    for h=1:h_horizon
        x_t_h = C'*(A'*B);  % gives coefficients of the Forecast X_(t+h) in terms of PCA basis
        if h < h_horizon
            C = (((Feigval(1:L)).^(-1)).* (x_t_h)')/(nobs-1);
        end
    end
    
    % calculating Prediction in terms of basis representation
    coef_eigf       = getcoef(PCA_Str.pcafd);
    coef_mean       = getcoef(PCA_Str.meanfd);
    FuncPredCoeff   = (x_t_h * coef_eigf(:,1:L)')'+coef_mean;

    
    
    
end
