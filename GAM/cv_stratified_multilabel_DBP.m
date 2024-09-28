function [test_indices train_indices ind]=cv_stratified_multilabel_DBP(N,n_folds,bestmodel)
% N is the label matrix
% Stratifies by discretizing multiclass tracking data into one unique label
% after binning individual classes. The final discretized label makes up a
% combination of individual classes labels.

% Bestmodel combinations:
% 1 = view, 2 = pl, 3 = hd, 4 = viewpl, 5 = viewhd, 6 = plhd, 7 = viewplhd
% xview 2 bins, yview 4 bins, zview 3 bins
% xplace 2 bins, yplace 4 bins, zplace 3 bins
% yaw 4 bins, pitch 3 bins, roll 2 bins

switch bestmodel

    case 1
        locx_view= discretize(N(:,1),2);          % Discretize x-view into 2 bins
        locy_view= discretize(N(:,2),4);          % Discretize y-view into 4 bins
        locz_view= discretize(N(:,3),3);          % Discretize z-view into 3 bins
        count=[length(unique(locx_view))  length(unique(locy_view))  length(unique(locz_view))];  % Get counts of unique bins for each axis
        ind = sub2ind(count,locx_view,locy_view,locz_view);  % Convert subscripts to linear indices

    case 2
        locx_place= discretize(N(:,1),2);          % Discretize x-place into 2 bins
        locy_place= discretize(N(:,2),4);          % Discretize y-place into 4 bins
        locz_place= discretize(N(:,3),3);          % Discretize z-place into 3 bins
        count=[length(unique(locx_place))  length(unique(locy_place))  length(unique(locz_place))];  % Get counts of unique bins for each axis
        ind = sub2ind(count,locx_place,locy_place,locz_place);  % Convert subscripts to linear indices
    case 3

        [counts, binEdges, locx_rot] = histcounts(N(:,1), [ -pi   -pi/2         0    pi/2    pi]);  % Discretize yaw into 4 bins
        [counts, binEdges, locy_rot] = histcounts(N(:,2), [-pi/2   -0.5236    0.5236    pi/2]);     % Discretize pitch into 3 bins
        [counts, binEdges, locz_rot] = histcounts(N(:,3),[-pi circ_mean(N(:,3)) pi]);               % Discretize roll into bins based on circular mean
        count=[length(unique(locx_rot))  length(unique(locy_rot))  length(unique(locz_rot))];        % Get counts of unique bins for each axis
        ind = sub2ind(count,locx_rot,locy_rot,locz_rot);  % Convert subscripts to linear indices

    case 4
        locx_view= discretize(N(:,1),2);          % Discretize x-view into 2 bins
        locy_view= discretize(N(:,2),4);          % Discretize y-view into 4 bins
        locz_view= discretize(N(:,3),3);          % Discretize z-view into 3 bins
        countgz=[length(unique(locx_view))  length(unique(locy_view))  length(unique(locz_view))];  % Counts for view axes
        locx_place= discretize(N(:,4),2);          % Discretize x-place into 2 bins
        locy_place= discretize(N(:,5),4);          % Discretize y-place into 4 bins
        locz_place= discretize(N(:,6),3);          % Discretize z-place into 3 bins
        countpl=[length(unique(locx_place))  length(unique(locy_place))  length(unique(locz_place))];  % Counts for place axes
        count=[countgz countpl];  % Combine counts
        ind = sub2ind(count,locx_view,locy_view,locz_view,locx_place,locy_place,locz_place);  % Convert subscripts to linear indices

    case 5
        locx_view= discretize(N(:,1),2);          % Discretize x-view into 2 bins
        locy_view= discretize(N(:,2),4);          % Discretize y-view into 4 bins
        locz_view= discretize(N(:,3),3);          % Discretize z-view into 3 bins
        countgz=[length(unique(locx_view))  length(unique(locy_view))  length(unique(locz_view))];  % Counts for view axes
        [counts, binEdges, locx_rot] = histcounts(N(:,4), [ -pi   -pi/2         0    pi/2    pi]);  % Discretize yaw into 4 bins
        [counts, binEdges, locy_rot] = histcounts(N(:,5), [-pi/2   -0.5236    0.5236    pi/2]);     % Discretize pitch into 3 bins
        [counts, binEdges, locz_rot] = histcounts(N(:,6),[-pi circ_mean(N(:,6)) pi]);               % Discretize roll into bins based on circular mean
        counthd=[length(unique(locx_rot))  length(unique(locy_rot))  length(unique(locz_rot))];      % Counts for rotation axes
        count=[countgz counthd];  % Combine counts
        ind = sub2ind(count,locx_view,locy_view,locz_view,locx_rot,locy_rot,locz_rot);  % Convert subscripts to linear indices

    case 6
        locx_place= discretize(N(:,1),2);          % Discretize x-place into 2 bins
        locy_place= discretize(N(:,2),4);          % Discretize y-place into 4 bins
        locz_place= discretize(N(:,3),3);          % Discretize z-place into 3 bins
        countpl=[length(unique(locx_place))  length(unique(locy_place))  length(unique(locz_place))];  % Counts for place axes
        [counts, binEdges, locx_rot] = histcounts(N(:,4), [ -pi   -pi/2         0    pi/2    pi]);     % Discretize yaw into 4 bins
        [counts, binEdges, locy_rot] = histcounts(N(:,5), [-pi/2   -0.5236    0.5236    pi/2]);        % Discretize pitch into 3 bins
        [counts, binEdges, locz_rot] = histcounts(N(:,6),[-pi circ_mean(N(:,6)) pi]);                  % Discretize roll into bins based on circular mean
        counthd=[length(unique(locx_rot))  length(unique(locy_rot))  length(unique(locz_rot))];         % Counts for rotation axes
        count=[countpl counthd];  % Combine counts
        ind = sub2ind(count,locx_place,locy_place,locz_place,locx_rot,locy_rot,locz_rot);  % Convert subscripts to linear indices

    case 7
        locx_view= discretize(N(:,1),2);          % Discretize x-view into 2 bins
        locy_view= discretize(N(:,2),4);          % Discretize y-view into 4 bins
        locz_view= discretize(N(:,3),3);          % Discretize z-view into 3 bins
        countgz=[length(unique(locx_view))  length(unique(locy_view))  length(unique(locz_view))];  % Counts for view axes
        locx_place= discretize(N(:,4),2);          % Discretize x-place into 2 bins
        locy_place= discretize(N(:,5),4);          % Discretize y-place into 4 bins
        locz_place= discretize(N(:,6),3);          % Discretize z-place into 3 bins
        countpl=[length(unique(locx_place))  length(unique(locy_place))  length(unique(locz_place))];  % Counts for place axes
        [counts, binEdges, locx_rot] = histcounts(N(:,7), [ -pi   -pi/2         0    pi/2    pi]);     % Discretize yaw into 4 bins
        [counts, binEdges, locy_rot] = histcounts(N(:,8), [-pi/2   -0.5236    0.5236    pi/2]);        % Discretize pitch into 3 bins
        [counts, binEdges, locz_rot] = histcounts(N(:,9),[-pi circ_mean(N(:,9)) pi]);                  % Discretize roll into bins based on circular mean
        counthd=[length(unique(locx_rot))  length(unique(locy_rot))  length(unique(locz_rot))];         % Counts for rotation axes
        count=[countgz countpl counthd];  % Combine counts
        ind = sub2ind(count,locx_view,locy_view,locz_view,locx_place,locy_place,locz_place,locx_rot,locy_rot,locz_rot);  % Convert subscripts to linear indices
end

cv = cvpartition(ind,"KFold",n_folds,"Stratify",true);  % Create stratified K-fold cross-validation partition

test_indices=zeros([length(N) n_folds]);  % Initialize test indices matrix
train_indices=zeros([length(N) n_folds]);  % Initialize training indices matrix
for i=1:n_folds
    test_indices(:,i)=cv.test(i);  % Assign test indices for fold i
    train_indices(:,i)=cv.training(i);  % Assign training indices for fold i
end

test_indices=logical(test_indices);  % Convert test indices to logical array
train_indices=logical(train_indices);  % Convert training indices to logical array
end
