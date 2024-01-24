%% This script runs the GAM analysis described in methods, and replicates figure S5c
% To install GAM_models, follow the instructions here (https://GAM_models.readthedocs.io/en/latest/) for pip installation
% agv models correspond to ANGULAR HEAD VELOCITY, bv to TRANSLATION SPEED 

clear all
close all
rng shuffle
load sample_data_int.mat
variable_types = repmat(cellstr('cell'), 1, 7);
shuffle_number=100;%how many shuffles
grid_number=100;%how many grid search lambda values
neuron='C20190530C21_12' ; %ID of Neuron to model
n_folds=5; % how many crossvalidation partitions
sr=60;% data sampling rate
jj=1;% neuron index
GAM_models=table;%stores the trained gam model
GAM_models.unitID{1}=neuron;
lams = rand([grid_number 1]);
lams = lams * 6 - 3 ;
lams = py.numpy.array([10.^lams]) ;
types_ID=string(compactWave.UnitIDs);



%Replace bellow with your python environment path
setenv('path',['...\AppData\Local\Programs\Python', getenv('path')]);
pe=pyenv('Version','...\AppData\Local\Programs\Python\Python310\pythonw.exe',...
    'ExecutionMode','InProcess');


%%
TrackingData=even_TData(TrackingData);
puls=TrackingData.neural_time;
fields=fieldnames(units.singleunits);
heatmrawf=Head_Orientation_No_floors(TrackingData);
invalid_frames=isnan(heatmrawf(:,1));
angl_vel=HeadSpeed(TrackingData,1/60,1);
body_vel=TracVelocity(TrackingData,0.0167,2);

%
unit=units.singleunits.(neuron);

unitID=neuron;
type_i=find(unitID==types_ID);
if isempty(type_i)
    error('no neuronal type info')

end
%%

neuron_type=compactWave.NeuronType(type_i);

gapt=find(gap)./30000; %in seconds

inv_frm=ismember(puls,gapt);
spikes=histc(unit,puls);
agv=angl_vel;
bdv=body_vel;

if length(spikes)> length(TrackingData.XPosition) %Cheks if Tracking data is longer than the neural data
    spikes=spikes(1:length(TrackingData.XPosition));
    agv(length(TrackingData.XPosition))=nan;
    bdv(length(TrackingData.XPosition))=nan;
    inv_frm=inv_frm(1:length(TrackingData.XPosition));
    invalid_frames=invalid_frames(1:length(TrackingData.XPosition));
else
    inv_frm=inv_frm(1:length(spikes),:);
    invalid_frames=invalid_frames(1:length(spikes));
    agv=agv(1:length(spikes),:);
    bdv=bdv(1:length(spikes),:);
end
invalid_frames(inv_frm)=true;

spikes(invalid_frames)=[];



agv(invalid_frames,:)=[];
bdv(invalid_frames,:)=[];


np=py.importlib.import_module('pygam');
agv(isnan(bdv))=nan;
spikes(isnan(agv))=[];
bdv(isnan(agv))=[];
agv(isnan(agv))=[];


[countt edges  loc] = histcounts(agv,20); %M x N MATRIX Histogram

%invalid Bins --> spent less than 1 second or visited less than 3
%times,assuming fps are 60
countt(countt<30)=nan;

[r c v]=size(countt);
binvisit=zeros(size(countt));
visit=abs(diff(loc));
index=find(visit);
binvisit(loc(index(1),1))=1;% first bin visited
binvisit(loc(index(end)+1))=1;% last bin visited
for i=1:length(find(visit))-1
    binvisit(loc(index(i+1),1))=binvisit(loc(index(i+1),1))+1;
end

countt(binvisit<3)=nan;

agv(ismember(loc,(find(isnan(countt)))))=nan;

[countt edges  loc] = histcounts(bdv,20); %M x N MATRIX Histogram

%invalid Bins --> spent less than 1 second or visited less than 3
%times,assuming fps are 60
countt(countt<30)=nan;

[r c v]=size(countt);
binvisit=zeros(size(countt));
visit=abs(diff(loc));
index=find(visit);
binvisit(loc(index(1),1))=1;% first bin visited
binvisit(loc(index(end)+1))=1;% last bin visited
for i=1:length(find(visit))-1
    binvisit(loc(index(i+1),1))=binvisit(loc(index(i+1),1))+1;
end

countt(binvisit<3)=nan;

bdv(ismember(loc,(find(isnan(countt)))))=nan;

spikes_agv=spikes;
spikes_bvd=spikes;

spikes_agv(isnan(agv))=[];

spikes_bvd(isnan(bdv))=[];
agv(isnan(bdv))=nan;

bdv(isnan(agv))=nan;
spikes_speed=spikes;
spikes_speed(isnan(bdv))=[];
agv(isnan(bdv))=[];
bdv(isnan(bdv))=[];


X_arrayang_vel=py.numpy.array([agv]);
X_arrayb_vel=py.numpy.array([bdv]);
X_arrayboth_vel=py.numpy.array([agv bdv]);

y_array_agv=py.numpy.array(spikes_speed);
y_array_bvd=py.numpy.array(spikes_speed);





if neuron_type=="Int"
    X_arrayang_vel=py.numpy.array([agv]);
    X_arrayb_vel=py.numpy.array([bdv]);
    X_arrayboth_vel=py.numpy.array([agv bdv]);
    gamagv=np.GAM(pyargs('distribution','poisson','link','log'));
    gamagv=gamagv.fit(X_arrayang_vel ,y_array_agv);
    gam_grid_agv=gamagv.gridsearch(X_arrayang_vel ,y_array_agv,pyargs('lam',lams));
    gam_grid_agv=pyGAM2struct(gam_grid_agv);
    lamagv=gam_grid_agv.lam;
    sr=60;
    shift=randi([5*sr length(spikes_speed)-(5*sr)],[shuffle_number,1]); %-+ 90 sec shift

    shuffle_agv=cell(shuffle_number,1);
    shuffle_ev_agv=nan(shuffle_number,1);
    for i=1:shuffle_number
        spikes_speed_sh=circshift(spikes_speed,shift(i));
        y_array_agv_sh=py.numpy.array(spikes_speed_sh);
        gam_grid_agv_shuffle=np.GAM(pyargs('distribution','poisson','link','log','lam',lamagv));
        gam_grid_agv_shuffle=gam_grid_agv_shuffle.fit(X_arrayang_vel ,y_array_agv_sh);
        gam_grid_agv_shuffle=pyGAM2struct(gam_grid_agv_shuffle);
        shuffle_agv{i}=gam_grid_agv_shuffle;
        shuffle_ev_agv(i)=gam_grid_agv_shuffle.pseudo_r2.explained_deviance;
    end


    % for i = 1:100
    %  shuffle_ev_agv(i)=   pyGAM.agv_shuffle{1,1}{i, 1}.pseudo_r2.explained_deviance;
    % end

    %statistics
    Rsh=shuffle_ev_agv;
    R=gam_grid_agv.pseudo_r2.explained_deviance;
    mu = mean(Rsh);
    sigma = std(Rsh);
    p_val = 1 - normcdf(R,mu,sigma);
    reject_h0 = p_val < 0.05/3;
    gam_grid_agv.p_val_shuffle=p_val;
    gam_grid_agv.reject_h0=reject_h0;
    gam_grid_agv.ED_normalized2shuffle=R-mu;
    gam_grid_agv.shuffle_ED_mu=mu;
    gam_grid_agv.shuffle_ED_sigma=sigma;


    gambv=np.GAM(pyargs('distribution','poisson','link','log'));
    gambv=gambv.fit(X_arrayb_vel,y_array_bvd);
    gam_grid_bv=gambv.gridsearch(X_arrayb_vel,y_array_bvd,pyargs('lam',lams));
    gam_grid_bv=pyGAM2struct(gam_grid_bv);
    lambv=gam_grid_bv.lam;

    shuffle_bv=cell(shuffle_number,1);
    shuffle_ev_bv=nan(shuffle_number,1);
    for i=1:shuffle_number
        spikes_speed_sh=circshift(spikes_speed,shift(i));
        y_array_bvd_sh=py.numpy.array(spikes_speed_sh);
        gam_grid_bv_sh=np.GAM(pyargs('distribution','poisson','link','log','lam',lambv));
        gam_grid_bv_sh=gam_grid_bv_sh.fit(X_arrayb_vel ,y_array_bvd_sh);
        gam_grid_bv_sh=pyGAM2struct(gam_grid_bv_sh);
        shuffle_bv{i}=gam_grid_bv_sh;
        shuffle_ev_bv(i)=gam_grid_bv_sh.pseudo_r2.explained_deviance;
    end

    Rsh=shuffle_ev_bv;
    R=gam_grid_bv.pseudo_r2.explained_deviance;
    mu = mean(Rsh);
    sigma = std(Rsh);
    p_val = 1 - normcdf(R,mu,sigma);
    reject_h0 = p_val < 0.05/3;
    gam_grid_bv.p_val_shuffle=p_val;
    gam_grid_bv.reject_h0=reject_h0;
    gam_grid_bv.ED_normalized2shuffle=R-mu;
    gam_grid_bv.shuffle_ED_mu=mu;
    gam_grid_bv.shuffle_ED_sigma=sigma;

    % for i = 1:100
    % shuffle_ev_bv(i)=   pyGAM.bv_shuffle{1,1}{i, 1}.pseudo_r2.explained_deviance;
    % end



    %

    models_reject=[gam_grid_agv.reject_h0 gam_grid_bv.reject_h0];
    models_ED=[gam_grid_agv.ED_normalized2shuffle gam_grid_bv.ED_normalized2shuffle];
    [m top_model]=max(models_ED)
    if gam_grid_agv.reject_h0 | gam_grid_bv.reject_h0

        if sum(models_reject)>1

            switch top_model

                case 1

                    [gam_grid_both, shuffle_both]=pyGAM_int_2models_shuffle(y_array_bvd,X_arrayboth_vel,2,np,[lamagv lambv],shift,shuffle_number)

                case 2
                    [gam_grid_both, shuffle_both]=pyGAM_int_2models_shuffle(y_array_bvd,X_arrayboth_vel,1,np,[lamagv lambv],shift,shuffle_number)

            end

        else

            x2shuffle=find(~models_reject);
            [gam_grid_both, shuffle_both]=pyGAM_int_2models_shuffle(y_array_bvd,X_arrayboth_vel,x2shuffle,np,[lamagv lambv],shift,shuffle_number);
        end
        GAM_models.agv_bv{jj}=gam_grid_both;
        GAM_models.agv_bv_shuffle{jj}=shuffle_both;
    end





    GAM_models.agv{jj}=gam_grid_agv;
    GAM_models.bv{jj}=gam_grid_bv;

    GAM_models.agv_shuffle{jj}=shuffle_agv;
    GAM_models.bv_shuffle{jj}=shuffle_bv;


else
    error('neuron is not a pyramidal cell')
end



%% Checking the best model
% bestmodel combinations 1=AHV 2=TS 3=AHV+TS
structGAM=GAM_models;
i=1; %neuron index
agv_dev_exp=[structGAM.agv{i,1}.ED_normalized2shuffle];
bv_dev_exp=[structGAM.bv{i,1}.ED_normalized2shuffle];
agv_reject=[structGAM.agv{i,1}.reject_h0];
bv_reject=structGAM.bv{i,1}.reject_h0;
if bv_reject|agv_reject %
    both_dev_exp=[structGAM.agv_bv{i,1}.ED_normalized2shuffle];
    both_dev_exp_real=[structGAM.agv_bv{i,1}.pseudo_r2.explained_deviance];
    agv_bv_mixed_reject=[structGAM.agv_bv{i,1}.reject_h0];
else
    agv_bv_mixed_reject=[nan];
    both_dev_exp=[nan];
end
bestmodelidx=[];
idint=[structGAM.unitID{i,1}];%neuron ID
mixed=find((agv_bv_mixed_reject>0));
bestmodelidx(mixed)=3;
ahvminusbv=agv_dev_exp-bv_dev_exp;
ahv_cells=find(ahvminusbv>0&(agv_bv_mixed_reject<1)&agv_reject);
bv_cells=find(ahvminusbv<0&(agv_bv_mixed_reject<1)&bv_reject);
bestmodelidx(ahv_cells)=1;
bestmodelidx(bv_cells)=2;




%% Stratified crossvalidation of the best model

bestmodel=bestmodelidx;
gamid=height(GAM_models);
j=height(GAM_models);
X_arrayang_vel=py.numpy.array([agv]);
X_arrayb_vel=py.numpy.array([bdv]);
X_arrayboth_vel=py.numpy.array([agv bdv]);
% lams = py.numpy.logspace(py.int(-3), py.int(3),py.int(3 ));

% [test_indices train_indices]=crossvalpartitionDBP(length(spikes_speed),n_folds);
spikes_model_total={};
spikes_real_total={};
modelcorr=nan(1,n_folds);
switch bestmodel
    case 1

    lams=GAM_models.agv{gamid,:}.lam;    
    [test_indices train_indices ind]=cv_stratified_multilabel_DBP_speed(agv,n_folds,bestmodel);  
    for nfold=1:n_folds
    y_array_agv=py.numpy.array(spikes_speed(train_indices(:,nfold)));
    X_arrayang_vel=py.numpy.array([agv(train_indices(:,nfold))]);    
    X_arrayang_vel_test=py.numpy.array([agv(test_indices(:,nfold))]);    
    spikes_test=spikes_speed(test_indices(:,nfold));
    spikes_real_total{nfold}=spikes_test;
    gamagv=np.GAM(pyargs('distribution','poisson','link','log','lam',lams));
    gamagv=gamagv.fit(X_arrayang_vel ,y_array_agv);
%   gam_grid_agv=gamagv.gridsearch(X_arrayang_vel ,y_array_agv);
    spikes_model=double(gamagv.predict(X_arrayang_vel_test));
    xr=corr(spikes_model',spikes_test);
    spikes_model_total{nfold}=spikes_model;
    modelcorr(nfold)=xr;
    end

    case 2
    lams=GAM_models.bv{gamid,:}.lam;    
    [test_indices train_indices ind]=cv_stratified_multilabel_DBP_speed(bdv,n_folds,bestmodel); 
    for nfold=1:n_folds
    y_array_agv=py.numpy.array(spikes_speed(train_indices(:,nfold)));
    X_arrayang_vel=py.numpy.array([bdv(train_indices(:,nfold))]);    
    X_arrayang_vel_test=py.numpy.array([bdv(test_indices(:,nfold))]);    
    spikes_test=spikes_speed(test_indices(:,nfold));
    spikes_real_total{nfold}=spikes_test;
    gamagv=np.GAM(pyargs('distribution','poisson','link','log','lam',lams));
    gamagv=gamagv.fit(X_arrayang_vel ,y_array_agv);
%   gam_grid_agv=gamagv.gridsearch(X_arrayang_vel ,y_array_agv);
    spikes_model=double(gamagv.predict(X_arrayang_vel_test));
    xr=corr(spikes_model',spikes_test);
    spikes_model_total{nfold}=spikes_model;
    modelcorr(nfold)=xr;
    end
   case 3
    lams=GAM_models.agv_bv{gamid,:}.lam;    
    [test_indices train_indices ind]=cv_stratified_multilabel_DBP_speed([agv bdv],n_folds,bestmodel); 
    for nfold=1:n_folds
    y_array_agv=py.numpy.array(spikes_speed(train_indices(:,nfold)));
    X_arrayang_vel=py.numpy.array([agv(train_indices(:,nfold)) bdv(train_indices(:,nfold))]);    
    X_arrayang_vel_test=py.numpy.array([agv(test_indices(:,nfold)) bdv(test_indices(:,nfold))]);    
    spikes_test=spikes_speed(test_indices(:,nfold));
    spikes_real_total{nfold}=spikes_test;
    gamagv=np.GAM(pyargs('distribution','poisson','link','log','lam',lams));
    gamagv=gamagv.fit(X_arrayang_vel ,y_array_agv);
    %gam_grid_agv=gamagv.gridsearch(X_arrayang_vel ,y_array_agv);
    spikes_model=double(gamagv.predict(X_arrayang_vel_test));
    xr=corr(spikes_model',spikes_test);
    spikes_model_total{nfold}=spikes_model;
    modelcorr(nfold)=xr;
    end
end


pyGAMcrossval.modelcorr{j}=modelcorr;
pyGAMcrossval.spikes_model_total{j}=spikes_model_total;
pyGAMcrossval.spikes_real_total{j}=spikes_real_total;
pyGAMcrossval.ind{j}=ind;
pyGAMcrossval.test_indices{j}=test_indices;
pyGAMcrossval.train_indices{j}=train_indices;


%% R squared F-test

p=[];
rejectint=[];
R2total=[];
int_selec_categ=string;
for i=1:height(GAM_models)
% % R2 F-test
ID=idint(i);
switch bestmodelidx(i)
    case 0
        int_selec_categ(i)="none";
        p(i)=nan;
        rejectint(i)=0;
        continue
    case 1

        R2=mean(pyGAMcrossval.modelcorr{j,1}).^2;
        k=floor(GAM_models.agv{j,1}.edof);
        int_selec_categ(i)="ahv";
    case 2

        R2=mean(pyGAMcrossval.modelcorr{j,1}).^2;
        k=floor(GAM_models.bv{j,1}.edof);
        int_selec_categ(i)="bv";
    case 3

        R2=mean(pyGAMcrossval.modelcorr{j,1}).^2;
        k=floor(GAM_models.agv_bv{j,1}.edof);
        int_selec_categ(i)="ahv+bv";
end
n=length(pyGAMcrossval.spikes_real_total{j,1}{1,1});
F = (R2/k) / ((1-R2)/(n-k-1));
alpha = 0.05;
Fcrit = finv(1-alpha,k,n-k-1);
p(i) = 1 - fcdf(F,k,n-k-1);
if p(i) < alpha
    rejectint(i)=1;
else
    rejectint(i)=0;
end
end
%% Plotting models predictions against firing rate
figure
X_arrayang_vel=py.numpy.array([agv bdv]);    
plot((1:length(spikes_speed))/sr,movmean(spikes_speed,sr/5).*sr,'Color',[0.4 0.4 0.4],'LineWidth',2);
hold on
spikes_model=double(gamagv.predict(X_arrayang_vel));
plot((1:length(spikes_speed))/sr,movmean(spikes_model,sr/5).*sr,'Color',[1 0.2 0.2],'LineWidth',1.5);
legend('real','model')
ylabel('Firing Rate Hz')
xlabel('Time (s)')
%% Plotting crossvalidation results

for i =1:n_folds
    figure
    plot(movmean(pyGAMcrossval.spikes_real_total{1, 1}{1, i},sr/5).*sr,'Color',[0.4 0.4 0.4],'LineWidth',2)
    hold on
    plot(movmean(pyGAMcrossval.spikes_model_total{1, 1}{1, i},sr/5).*sr,'Color',[1 0.2 0.2],'LineWidth',1.5)
legend('real','model')
end




