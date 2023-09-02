%%
clear all
close all
rng shuffle
load sample_data.mat


pyGAM=table;
types_ID=string(compactWave.UnitIDs);

   
TrackingData=even_TData(TrackingData);


puls=TrackingData.neural_time;
fields=fieldnames(units.singleunits);

heatmrawf=Head_Orientation_No_floors(TrackingData);


invalid_frames=isnan(heatmrawf(:,1));
angl_vel=HeadSpeed(TrackingData,1/60,1);

body_vel=TracVelocity(TrackingData,0.0167,2);

temporal=table;
for j=1:length(fields)
%%
neuron=char(fields(j))

unit=units.singleunits.(neuron);
    
unitID=string(fields(j));
type_i=find(unitID==types_ID);
if isempty(type_i) 
    continue
end
%%
neuron_type=compactWave.NeuronType(type_i);

temporal.unitID{j}=unitID;


%removing first optirack frames because they are not recorded by optitrack
%finding number of spikes per frame


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

% spikes=movmean(spikes,3);


% yedge=[0 0.2 0.4 0.6 0.8 1 1.2];
% xedge=[0 0.2 0.4 0.61];
% zedge=[ 0 0.2 0.4 0.6 0.8 1 1.2 1.41];%HEIGHT
% 
% [count edges mid loc] = histcn([Xgf Ygf Zgf],5,3,5); %M x N MATRIX Histogram
% 
% gaze_cat=sub2ind(size(count),loc(:,1),loc(:,2),loc(:,3));
% 
% [count edges mid loc] = histcn([X Y Z],5,3,5); %M x N MATRIX Histogram
% 
% place_cat=sub2ind(size(count),loc(:,1),loc(:,2),loc(:,3));

close all

%% make sure you install scipy and numpy with pip , not the dfault anaconda libraries
% setenv('path',['C:\Users\diego\anaconda3', getenv('path')]);
% pe=pyenv('Version',['C:\Users\diego\anaconda3\pythonw.exe'],...
%      'ExecutionMode','InProcess');

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

% cv=struct;
% test_size=ceil(length(agv).*0.3);
% test=randi([1 length(agv)-test_size]);
% test=test:1:test+test_size;
% cv.test=false([length(agv) 1]);
% cv.train=false([length(agv) 1]);
% cv.test(test)=true;
% cv.train(~cv.test)=true;

X_arrayang_vel=py.numpy.array([agv]);
X_arrayb_vel=py.numpy.array([bdv]);
X_arrayboth_vel=py.numpy.array([agv bdv]);

y_array_agv=py.numpy.array(spikes_speed);
y_array_bvd=py.numpy.array(spikes_speed);





% null_array=py.numpy.array(ones(size(double(X_array))));dtype=['numerical', 'categorical', 'numerical'
if neuron_type=="Int"

% cv=struct;
% test_size=ceil(length(agv).*0.2);
% test= 1%test_size*1;%randi([1 length(agv)-test_size]);
% 
% 
% test=test:1:test+test_size;
% cv.test=false([length(agv) 1]);
% cv.train=false([length(agv) 1]);
% cv.test(test)=true;
% cv.train(~cv.test)=true;


X_arrayang_vel=py.numpy.array([agv]);
X_arrayb_vel=py.numpy.array([bdv]);
X_arrayboth_vel=py.numpy.array([agv bdv]);
lams = py.numpy.logspace(py.int(-3), py.int(3),py.int(3 ));


gamagv=np.GAM(pyargs('distribution','poisson','link','log'));
gamagv=gamagv.fit(X_arrayang_vel ,y_array_agv);
% gam_grid_agv=gamagv.gridsearch(X_arrayang_vel ,y_array_agv,pyargs('lam',lams));
gam_grid_agv=pyGAM2struct(gamagv);
sr=60;
shift=randi([5*sr length(spikes_speed)-(5*60)],[50,1]); %-+ 90 sec shift

shuffle_agv=cell(50,1);
shuffle_ev_agv=nan(50,1);
for i=1:50
spikes_speed_sh=circshift(spikes_speed,shift(i)); 
y_array_agv_sh=py.numpy.array(spikes_speed_sh);
gamagv=np.GAM(pyargs('distribution','poisson','link','log'));
gamagv=gamagv.fit(X_arrayang_vel ,y_array_agv_sh);
% gam_grid_agv_shuffle=gamagv.gridsearch(X_arrayang_vel ,y_array_agv_sh,pyargs('lam',lams));
gam_grid_agv_shuffle=pyGAM2struct(gamagv);
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
gambv=gambv.fit(X_arrayb_vel ,y_array_bvd);
% gam_grid_bv=gambv.gridsearch(X_arrayb_vel,y_array_bvd,pyargs('lam',lams));
gam_grid_bv=pyGAM2struct(gambv);

shuffle_bv=cell(50,1);
shuffle_ev_bv=nan(50,1);
for i=1:50
spikes_speed_sh=circshift(spikes_speed,shift(i)); 
y_array_bvd_sh=py.numpy.array(spikes_speed_sh);
gambv=np.GAM(pyargs('distribution','poisson','link','log'));
gambv=gambv.fit(X_arrayb_vel ,y_array_bvd_sh);
% gam_grid_bv_sh=gambv.gridsearch(X_arrayb_vel,y_array_bvd_sh,pyargs('lam',lams));
gam_grid_bv_sh=pyGAM2struct(gambv);
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

             [gam_grid_both, shuffle_both]=pyGAM_int_2models_shuffle(y_array_bvd,X_arrayboth_vel,2,np,lams,shift);

           case 2
             [gam_grid_both, shuffle_both]=pyGAM_int_2models_shuffle(y_array_bvd,X_arrayboth_vel,1,np,lams,shift);

       end

    else
    
    x2shuffle=find(~models_reject);    
    [gam_grid_both, shuffle_both]=pyGAM_int_2models_shuffle(y_array_bvd,X_arrayboth_vel,x2shuffle,np,lams,shift);
    end
    temporal.agv_bv{j}=gam_grid_both;
    temporal.agv_bv_shuffle{j}=shuffle_both;
end





temporal.agv{j}=gam_grid_agv;
temporal.bv{j}=gam_grid_bv;

temporal.agv_shuffle{j}=shuffle_agv;
temporal.bv_shuffle{j}=shuffle_bv;


else
    continue
end

end

if any(strcmp('agv',fieldnames(temporal)))
pyGAM=[pyGAM;temporal];
end




