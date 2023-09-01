function [gam_grid_both, shuffle_both] = pyGAM_pyr_nogrid_shuffle(y_array_bvd, X_arrayboth_vel, x2shuffle, np, lams, shift)
    % Perform Poisson GAM fitting and shuffling analysis
    
    % Create a GAM model for original data
    gamboth = np.GAM(pyargs('distribution', 'poisson', 'link', 'log'));
    gamboth = gamboth.fit(X_arrayboth_vel, y_array_bvd);
    
    % Convert the GAM model to a structured grid
    gam_grid_both = pyGAM2struct(gamboth);
    
    % Initialize storage for shuffled results
    shuffle_both = cell(50, 1);
    shuffle_ED_both = nan(50, 1);
    
    % Perform shuffling analysis
    for i = 1:50
        x_arr = double(X_arrayboth_vel);
        x_arr(:, x2shuffle) = circshift(x_arr(:, x2shuffle), shift(i));
        x_array_both_sh = py.numpy.array(x_arr);
        
        % Create a GAM model for shuffled data
        gamboth_sh = np.GAM(pyargs('distribution', 'poisson', 'link', 'log'));
        gamboth_sh = gamboth_sh.fit(x_array_both_sh, y_array_bvd);
        
        % Convert the shuffled GAM model to a structured grid
        gam_grid_both_sh = pyGAM2struct(gamboth_sh);
        
        % Store shuffled results and explained deviance
        shuffle_both{i} = gam_grid_both_sh;
        shuffle_ED_both(i) = gam_grid_both_sh.pseudo_r2.explained_deviance;
    end
    
    % Calculate statistics for shuffled results
    Rsh = shuffle_ED_both;
    R = gam_grid_both.pseudo_r2.explained_deviance;
    mu = mean(Rsh);
    sigma = std(Rsh);
    p_val = 1 - normcdf(R, mu, sigma);
    reject_h0 = p_val < 0.05 / 3;
    
    % Store statistical results in the structured grid
    gam_grid_both.p_val_shuffle = p_val;
    gam_grid_both.reject_h0 = reject_h0;
    gam_grid_both.ED_normalized2shuffle = R - mu;
    gam_grid_both.shuffle_ED_mu = mu;
    gam_grid_both.shuffle_ED_sigma = sigma;
    gam_grid_both.predictor2shuffle = x2shuffle;
end
