function [test_indices train_indices ind]=cv_stratified_multilabel_DBP_speed(N,n_folds,bestmodel)
% N is the label matrix
%Stratifies by discretizing multiclass tracking data into one unique label
%after bining individual classes. The final discretized label makes up an
%combination of individual classes labels.


if bestmodel==1||bestmodel==2
        ind= discretize(N,4);          
        
elseif bestmodel==3
        [count edges mid loc] =histcn(N,4,4);
        ind = sub2ind(size(count),loc(:,1),loc(:,2));
end

cv = cvpartition(ind,"KFold",n_folds,"Stratify",true);

test_indices=zeros([length(N) n_folds]);
train_indices=zeros([length(N) n_folds]);
for i=1:n_folds
test_indices(:,i)=cv.test(i);
train_indices(:,i)=cv.training(i);
end

test_indices=logical(test_indices);
train_indices=logical(train_indices);
end


