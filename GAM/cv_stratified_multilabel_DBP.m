function [test_indices train_indices ind]=cv_stratified_multilabel_DBP(N,n_folds,bestmodel)
% N is the label matrix
%Stratifies by discretizing multiclass tracking data into one unique label
%after bining individual classes. The final discretized label makes up an
%combination of individual classes labels.

%bestmodel combinations 1=view 2=pl 3=hd 4=viewpl 5=viewhd 6= plhd 7= viewplhd
%xview 2 bins, y view 4 bins, zview 3 bins
%xplace 2 bins, y place 4 bins, zplace 3 bins
%yaw 4 bins pitch 3 bins roll 2 bins
switch bestmodel

    case 1
        locx_view= discretize(N(:,1),2);          
        locy_view= discretize(N(:,2),4); 
        locz_view= discretize(N(:,3),3); 
        count=[length(unique(locx_view))  length(unique(locy_view))  length(unique(locz_view))];
        ind = sub2ind(count,locx_view,locy_view,locz_view);

    case 2
        locx_place= discretize(N(:,1),2);          
        locy_place= discretize(N(:,2),4); 
        locz_place= discretize(N(:,3),3); 
        count=[length(unique(locx_place))  length(unique(locy_place))  length(unique(locz_place))];
        ind = sub2ind(count,locx_place,locy_place,locz_place);
    case 3
    
        [counts, binEdges, locx_rot] = histcounts(N(:,1), [ -pi   -pi/2         0    pi/2    pi]);
        [counts, binEdges, locy_rot] = histcounts(N(:,2), [-pi/2   -0.5236    0.5236    pi/2]);
        [counts, binEdges, locz_rot] = histcounts(N(:,3),[-pi circ_mean(N(:,3)) pi]); 
        count=[length(unique(locx_rot))  length(unique(locy_rot))  length(unique(locz_rot))];
        ind = sub2ind(count,locx_rot,locy_rot,locz_rot);

    case 4
        locx_view= discretize(N(:,1),2);          
        locy_view= discretize(N(:,2),4); 
        locz_view= discretize(N(:,3),3); 
        countgz=[length(unique(locx_view))  length(unique(locy_view))  length(unique(locz_view))];
        locx_place= discretize(N(:,4),2);          
        locy_place= discretize(N(:,5),4); 
        locz_place= discretize(N(:,6),3); 
        countpl=[length(unique(locx_place))  length(unique(locy_place))  length(unique(locz_place))];
        count=[countgz countpl];
        ind = sub2ind(count,locx_view,locy_view,locz_view,locx_place,locy_place,locz_place);
    
    case 5
        locx_view= discretize(N(:,1),2);          
        locy_view= discretize(N(:,2),4); 
        locz_view= discretize(N(:,3),3); 
        countgz=[length(unique(locx_view))  length(unique(locy_view))  length(unique(locz_view))];
        [counts, binEdges, locx_rot] = histcounts(N(:,4), [ -pi   -pi/2         0    pi/2    pi]);
        [counts, binEdges, locy_rot] = histcounts(N(:,5), [-pi/2   -0.5236    0.5236    pi/2]);
        [counts, binEdges, locz_rot] = histcounts(N(:,6),[-pi circ_mean(N(:,6)) pi]); 
        counthd=[length(unique(locx_rot))  length(unique(locy_rot))  length(unique(locz_rot))];
        count=[countgz counthd];
        ind = sub2ind(count,locx_view,locy_view,locz_view,locx_rot,locy_rot,locz_rot);
    
    case 6
        locx_place= discretize(N(:,1),2);          
        locy_place= discretize(N(:,2),4); 
        locz_place= discretize(N(:,3),3); 
        countpl=[length(unique(locx_place))  length(unique(locy_place))  length(unique(locz_place))];
        [counts, binEdges, locx_rot] = histcounts(N(:,4), [ -pi   -pi/2         0    pi/2    pi]);
        [counts, binEdges, locy_rot] = histcounts(N(:,5), [-pi/2   -0.5236    0.5236    pi/2]);
        [counts, binEdges, locz_rot] = histcounts(N(:,6),[-pi circ_mean(N(:,6)) pi]); 
        counthd=[length(unique(locx_rot))  length(unique(locy_rot))  length(unique(locz_rot))];
        count=[countpl counthd];
        ind = sub2ind(count,locx_place,locy_place,locz_place,locx_rot,locy_rot,locz_rot);
    
    case 7
        locx_view= discretize(N(:,1),2);          
        locy_view= discretize(N(:,2),4); 
        locz_view= discretize(N(:,3),3); 
        countgz=[length(unique(locx_view))  length(unique(locy_view))  length(unique(locz_view))];
        locx_place= discretize(N(:,4),2);          
        locy_place= discretize(N(:,5),4); 
        locz_place= discretize(N(:,6),3); 
        countpl=[length(unique(locx_place))  length(unique(locy_place))  length(unique(locz_place))];
        [counts, binEdges, locx_rot] = histcounts(N(:,7), [ -pi   -pi/2         0    pi/2    pi]);
        [counts, binEdges, locy_rot] = histcounts(N(:,8), [-pi/2   -0.5236    0.5236    pi/2]);
        [counts, binEdges, locz_rot] = histcounts(N(:,9),[-pi circ_mean(N(:,9)) pi]); 
        counthd=[length(unique(locx_rot))  length(unique(locy_rot))  length(unique(locz_rot))];
        count=[countgz countpl counthd];
        ind = sub2ind(count,locx_view,locy_view,locz_view,locx_place,locy_place,locz_place,locx_rot,locy_rot,locz_rot);
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


