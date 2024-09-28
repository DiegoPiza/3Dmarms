function [test_indices train_indices ind]=cv_stratified_multilabel_DBP_speed(N,n_folds,bestmodel)
% N is the label matrix
% Stratifies by discretizing multiclass tracking data into one unique label
% after binning individual classes. The final discretized label makes up a
% combination of individual classes labels.

if bestmodel==1||bestmodel==2      % Check if bestmodel is 1 or 2
        ind= discretize(N,4);          % Discretize N into 4 bins and assign to ind

elseif bestmodel==3                % If bestmodel is 3
        [count edges mid loc] =histcn(N,4,4);    % Compute multidimensional histogram counts
        ind = sub2ind(size(count),loc(:,1),loc(:,2));  % Convert subscripts to linear indices
end

cv = cvpartition(ind,"KFold",n_folds,"Stratify",true);    % Create a stratified K-fold partition

test_indices=zeros([length(N) n_folds]);       % Initialize test_indices matrix
train_indices=zeros([length(N) n_folds]);      % Initialize train_indices matrix
for i=1:n_folds                                % Loop over each fold
test_indices(:,i)=cv.test(i);                  % Get test indices for fold i
train_indices(:,i)=cv.training(i);             % Get training indices for fold i
end

test_indices=logical(test_indices);            % Convert test_indices to logical array
train_indices=logical(train_indices);          % Convert train_indices to logical array
end
