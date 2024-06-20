
%% SVM Place decoding
% Load the sample data from the provided MAT-file
load ensemble_sample_data.mat

% Initialize a matrix to store accuracy results, with dimensions corresponding to the third and second dimensions of the input data
accuracytotal=nan([size(all_sessionsfloor_iteration,3) size(all_sessionsfloor_iteration,2)]);

% Loop through each session (third dimension)
for ii=1:size(all_sessionsfloor_iteration,3)
    parfor jj=1:size(all_sessionsfloor_iteration,2)    
    [trainedClassifier, validationAccuracy] = trainClassifierfloor(all_sessionsfloor_iteration(:,jj,ii), modifiedfloor);
    accuracytotal(ii,jj)=validationAccuracy;
    end
end
%select best model
accuracytotal=mean(accuracytotal);
[best_decoding best_index]=max(accuracytotal);
best_decoding_total(1)=best_decoding;
best_index_total(1)=best_index;
selectedFeatures(:,:,:)=all_sessionsfloor_iteration(:,best_index,:);

% Loop through ensemble sizes from 2 to 20
for ens_size=2:20
% Create an index array excluding the best index from the previous iteration
    jjindex=1:size(all_sessionsfloor_iteration,2);
    jjindex(best_index_total)=[]; 
    accuracytotaltemp=nan([size(all_sessionsfloor_iteration,3) size(jjindex,2)]);
    accuracytotal=nan([size(all_sessionsfloor_iteration,3) size(all_sessionsfloor_iteration,2)]);
    for ii=1:size(all_sessionsfloor_iteration,3)
       
        parfor jj=1:length(jjindex)
        features=[selectedFeatures(:,:,ii) all_sessionsfloor_iteration(:,jjindex(jj),ii)];    
        [trainedClassifier, validationAccuracy] = trainClassifierfloor(features, modifiedfloor);
        %cellindex=jjindex(jj);
        accuracytotaltemp(ii,jj)=validationAccuracy;
        end
    end
% Copy the temporary accuracy results to the main accuracy matrix
accuracytotal(:,jjindex)=accuracytotaltemp;
% Compute the mean accuracy across all sessions
accuracytotal=nanmean(accuracytotal);
% Identify the best model based on the highest accuracy
[best_decoding best_index]=max(accuracytotal);
% Store the best decoding accuracy and corresponding index for the current ensemble size
best_decoding_total(ens_size)=best_decoding;
% Update the selected features with the best feature set found in the current iteration
best_index_total(ens_size)=best_index;
selectedFeatures=[selectedFeatures all_sessionsfloor_iteration(:,best_index,:)];
end
