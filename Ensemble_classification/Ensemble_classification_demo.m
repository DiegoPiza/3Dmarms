
%% Place decoding

load ensemble_sample_data.mat

accuracytotal=nan([size(all_sessionsfloor_iteration,3) size(all_sessionsfloor_iteration,2)]);

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
%
for ens_size=2:20
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
accuracytotal(:,jjindex)=accuracytotaltemp;
accuracytotal=nanmean(accuracytotal);
%select best model
[best_decoding best_index]=max(accuracytotal);
best_decoding_total(ens_size)=best_decoding;
best_index_total(ens_size)=best_index;
selectedFeatures=[selectedFeatures all_sessionsfloor_iteration(:,best_index,:)];
save('best20unitsC_50trials_4place1346Bin','selectedFeatures','best_index_total','best_decoding_total','neuronsIDTotalfl')
end
