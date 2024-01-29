function [gam_grid_both, shuffle_both]=pyGAM_int_2models_shuffle(y_array_bvd,X_arrayboth_vel,x2shuffle,np,lams,shift,shuffle_number)
    % pyGAM_int_2models_shuffle performs GAM analysis with shuffling for interaction models.
    % This function fits GAM models, shuffles specific predictor columns,
    % and evaluates the effect of shuffling on model performance.
    %
    % Inputs:
    %   y_array_bvd: Response variable data
    %   X_arrayboth_vel: Predictor variable data (both predictors)
    %   x2shuffle: Predictor columns to shuffle
    %   np: Python numpy module
    %   lams: Lambda values for the grid search
    %   shift: Array of shift values for shuffling
    %
    % Outputs:
    %   gam_grid_both: Struct containing GAM model grid search results
    %   shuffle_both: Cell array containing shuffled GAM model results
    
    % Fit the initial GAM model to the unshuffled data
gamboth=np.GAM(pyargs('distribution','poisson','link','log','lam',lams));
gamboth=gamboth.fit(X_arrayboth_vel ,y_array_bvd);
    % Initialize arrays for shuffled results

gam_grid_both=pyGAM2struct(gamboth);
shuffle_both=cell(shuffle_number,1);
shuffle_ED_both=nan(shuffle_number,1);
    % Loop for shuffling and fitting shuffled GAM models

for i=1:shuffle_number
x_arr=double(X_arrayboth_vel);
x_arr(:,x2shuffle)=circshift(x_arr(:,x2shuffle),shift(i)); 
x_array_both_sh=py.numpy.array(x_arr);
gamboth=np.GAM(pyargs('distribution','poisson','link','log','lam',lams));
gamboth=gamboth.fit(x_array_both_sh,y_array_bvd);
gam_grid_both_sh=pyGAM2struct(gamboth);
shuffle_both{i}=gam_grid_both_sh;
shuffle_ED_both(i)=gam_grid_both_sh.pseudo_r2.explained_deviance;
end

    % Calculate p-value and test hypothesis of shuffled effect

    Rsh=shuffle_ED_both;
    R=gam_grid_both.pseudo_r2.explained_deviance;
    mu = mean(Rsh);
        if isnan(mu)
        error('nan mean detected')
    end
    sigma = std(Rsh);
    p_val = 1 - normcdf(R,mu,sigma);
    reject_h0 = p_val < 0.05/3;
    gam_grid_both.p_val_shuffle=p_val;
    gam_grid_both.reject_h0=reject_h0;
    gam_grid_both.ED_normalized2shuffle=R-mu;
    gam_grid_both.shuffle_ED_mu=mu;
    gam_grid_both.shuffle_ED_sigma=sigma;
    gam_grid_both.predictor2shuffle=x2shuffle;
end
