function structGAM = pyGAM2struct(pGAM)

    % Extract statistics from pyGAM object to save as MATLAB file
    stats = struct(pGAM.statistics_);
    stats.pseudo_r2 = struct(stats.pseudo_r2);
    stats.GCV = [];
    stats.cov = [];
    stats.edof_per_coef = []; % struct(stats.edof_per_coef);
    stats.m_features = int64(stats.m_features);
    stats.n_samples = int64(stats.n_samples);
    stats.se = double(stats.se);
    stats.p_values = double(py.numpy.array(stats.p_values));
    stats.coef_ = double(pGAM.coef_);
    
    % Store statistics in structGAM
    structGAM = stats;
end
