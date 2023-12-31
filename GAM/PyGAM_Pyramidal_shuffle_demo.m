%% This script runs the GAM analysis described in methods, figure S3a, on putative pyramidal neurons
% To install PyGAM, follow the instructions here (https://pygam.readthedocs.io/en/latest/) for pip installation 
clear all
close all
rng shuffle
load sample_data.mat

types_ID=string(compactWave.UnitIDs);
%% 


TrackingData=even_TData(TrackingData);

puls=TrackingData.neural_time;
fields=fieldnames(units.singleunits);

heatmrawf=Head_Orientation_No_floors(TrackingData);

heatm=[];
heatm(:,1)=TrackingData.XPosition;
heatm(:,2)=TrackingData.YPosition;
heatm(:,3)=TrackingData.ZPosition;
invalid_frames=isnan(heatmrawf(:,1));
heatmpos=heatm;

q1=TrackingData.WQuat;
q2=TrackingData.iXQuat;
q3=TrackingData.jYQuat;
q4=TrackingData.kZQuat;

[x_rot y_rot z_rot] = quat2angle([q1 q2 q3 q4]);

temporal=table('Size',[1,15],'VariableTypes',variable_types,'VariableNames',variable_names);

for j=1:length(fields)
neuron=char(fields(j))

unit=units.singleunits.(neuron);
    
unitID=string(fields(j));
type_i=find(unitID==types_ID);
if isempty(type_i) 
    continue
end
neuron_type=compactWave.NeuronType(type_i);

temporal.unitID{j}=unitID;


%removing first optirack frames because they are not recorded by optitrack
%finding number of spikes per frame


gapt=find(gap)./30000; %in seconds

inv_frm=ismember(puls,gapt);
spikes=histc(unit,puls);
heatm=heatmpos;     
heatmgzf=heatmrawf;


if length(spikes)> length(TrackingData.XPosition) %Cheks if Tracking data is longer than the neural data
spikes=spikes(1:length(TrackingData.XPosition));
inv_frm=inv_frm(1:length(TrackingData.XPosition));
invalid_frames=invalid_frames(1:length(TrackingData.XPosition));
else 
heatm=heatm(1:length(spikes),:); 
inv_frm=inv_frm(1:length(spikes),:);
invalid_frames=invalid_frames(1:length(spikes));
heatmgzf=heatmgzf(1:length(spikes),:);
end
invalid_frames(inv_frm)=true;

spikes(invalid_frames)=[];
heatm(invalid_frames,:)=[];

%%%Fitting glm to position data


X=heatm(:,1);
Y=heatm(:,2);
Z=heatm(:,3);




heatmgzf(invalid_frames,:)=[];
Xgf=heatmgzf(:,1);
Ygf=heatmgzf(:,2);
Zgf=heatmgzf(:,3);

x_r = x_rot;
y_r = y_rot;
z_r = z_rot;

x_r(invalid_frames,:)=[];
y_r(invalid_frames,:)=[];
z_r(invalid_frames,:)=[];



% spikes=movmean(spikes,3);
% Xgf=movmean(Xgf,3);         
% Ygf=movmean(Ygf,3);
% Zgf=movmean(Zgf,3);
% 
% X=movmean(X,3);
% Y=movmean(Y,3);
% Z=movmean(Z,3);


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
% setenv('path',['C:\Users\dbuitra2\AppData\Local\Programs\Python\Python37', getenv('path')]);
% pe=pyenv('Version',['C:\Users\dbuitra2\AppData\Local\Programs\Python\Python37\pythonw.exe'],...
%      'ExecutionMode','InProcess');

np=py.importlib.import_module('pygam');
spikes(isnan(Zgf))=[];
X(isnan(Zgf))=[];
Y(isnan(Zgf))=[];
Z(isnan(Zgf))=[];
x_r(isnan(Zgf))=[];
y_r(isnan(Zgf))=[];
z_r(isnan(Zgf))=[];
Xgf(isnan(Zgf))=[];
Ygf(isnan(Zgf))=[];
Zgf(isnan(Zgf))=[];



% cv=struct;
% test_size=ceil(length(agv).*0.3);
% test=randi([1 length(agv)-test_size]);
% test=test:1:test+test_size;
% cv.test=false([length(agv) 1]);
% cv.train=false([length(agv) 1]);
% cv.test(test)=true;
% cv.train(~cv.test)=true;
X_arrayhd=py.numpy.array([x_r y_r z_r]);
X_arraygz=py.numpy.array([Xgf Ygf Zgf]);
X_arrayplace=py.numpy.array([X Y Z]);
 
X_arraygz_place=py.numpy.array([Xgf Ygf Zgf X Y Z]);

X_array_gz_place_hd=py.numpy.array([Xgf Ygf Zgf X Y Z x_r y_r z_r]);

X_array_gz_hd=py.numpy.array([Xgf Ygf Zgf x_r y_r z_r]);

X_array_place_hd=py.numpy.array([X Y Z x_r y_r z_r]);

X_array_gz_hd=py.numpy.array([Xgf Ygf Zgf x_r y_r z_r]);


y_array=py.numpy.array(spikes);


% null_array=py.numpy.array(ones(size(double(X_array))));dtype=['numerical', 'categorical', 'numerical'

if neuron_type=="Pyr"


lams = py.numpy.logspace(py.int(-3), py.int(3),py.int(3 ));
gamgz=np.GAM(pyargs('distribution','poisson','link','log'));
gamgz=gamgz.fit(X_arraygz ,y_array);
% gam_grid_gz=gamgz.gridsearch(X_arraygz ,y_array,pyargs('lam',lams));
gam_grid_gz=pyGAM2struct(gamgz);

sr=20;
shift=randi([5*sr length(spikes)-(5*sr)],[50,1]); %-+ 90 sec shift

shuffle_gz=cell(50,1);
shuffle_ev_gz=nan(50,1);
for i=1:50
spikes_sh=circshift(spikes,shift(i)); 
y_array_gz_sh=py.numpy.array(spikes_sh);
gamgz_sh=np.GAM(pyargs('distribution','poisson','link','log'));
gamgz_sh=gamgz_sh.fit(X_arraygz ,y_array_gz_sh);
% gam_grid_gz_shuffle=gamgz.gridsearch(X_arraygz ,y_array_gz_sh,pyargs('lam',lams));
gam_grid_gz_shuffle=pyGAM2struct(gamgz_sh);
shuffle_gz{i}=gam_grid_gz_shuffle;
shuffle_ev_gz(i)=gam_grid_gz_shuffle.pseudo_r2.explained_deviance;
end

%statistics
    Rsh=shuffle_ev_gz;
    R=gam_grid_gz.pseudo_r2.explained_deviance;
    mu = mean(Rsh);
    sigma = std(Rsh);
    p_val = 1 - normcdf(R,mu,sigma);
    reject_h0 = p_val < 0.05/6;
    gam_grid_gz.p_val_shuffle=p_val;
    gam_grid_gz.reject_h0=reject_h0;
    gam_grid_gz.ED_normalized2shuffle=R-mu;
    gam_grid_gz.shuffle_ED_mu=mu;
    gam_grid_gz.shuffle_ED_sigma=sigma;



gampl=np.GAM(pyargs('distribution','poisson','link','log'));
gampl=gampl.fit(X_arrayplace ,y_array);
% gam_grid_pl=gampl.gridsearch(X_arrayplace ,y_array,pyargs('lam',lams));
gam_grid_pl=pyGAM2struct(gampl);


shuffle_pl=cell(50,1);
shuffle_ev_pl=nan(50,1);
for i=1:50
spikes_sh=circshift(spikes,shift(i)); 
y_array_pl_sh=py.numpy.array(spikes_sh);
gampl_sh=np.GAM(pyargs('distribution','poisson','link','log'));
gampl_sh=gampl_sh.fit(X_arrayplace ,y_array_pl_sh);
% gam_grid_pl_shuffle=gampl.gridsearch(X_arrayplace ,y_array_pl_sh,pyargs('lam',lams));
gam_grid_pl_shuffle=pyGAM2struct(gampl_sh);
shuffle_pl{i}=gam_grid_pl_shuffle;
shuffle_ev_pl(i)=gam_grid_pl_shuffle.pseudo_r2.explained_deviance;
end

%statistics
    Rsh=shuffle_ev_pl;
    R=gam_grid_pl.pseudo_r2.explained_deviance;
    mu = mean(Rsh);
    sigma = std(Rsh);
    p_val = 1 - normcdf(R,mu,sigma);
    reject_h0 = p_val < 0.05/6;
    gam_grid_pl.p_val_shuffle=p_val;
    gam_grid_pl.reject_h0=reject_h0;
    gam_grid_pl.ED_normalized2shuffle=R-mu;
    gam_grid_pl.shuffle_ED_mu=mu;
    gam_grid_pl.shuffle_ED_sigma=sigma;


gamhd=np.GAM(pyargs('distribution','poisson','link','log'));
gamhd=gamhd.fit(X_arrayhd ,y_array);
% gam_grid_hd=gamhd.gridsearch(X_arrayhd ,y_array,pyargs('lam',lams));
gam_grid_hd=pyGAM2struct(gamhd);


shuffle_hd=cell(50,1);
shuffle_ev_hd=nan(50,1);
for i=1:50
spikes_sh=circshift(spikes,shift(i)); 
y_array_hd_sh=py.numpy.array(spikes_sh);
gamhd_sh=np.GAM(pyargs('distribution','poisson','link','log'));
gamhd_sh=gamhd_sh.fit(X_arrayhd ,y_array_hd_sh);
% gam_grid_hd_shuffle=gamhd.gridsearch(X_arrayhd ,y_array_hd_sh,pyargs('lam',lams));
gam_grid_hd_shuffle=pyGAM2struct(gamhd_sh);
shuffle_hd{i}=gam_grid_hd_shuffle;
shuffle_ev_hd(i)=gam_grid_hd_shuffle.pseudo_r2.explained_deviance;
end

%statistics
    Rsh=shuffle_ev_hd;
    R=gam_grid_hd.pseudo_r2.explained_deviance;
    mu = mean(Rsh);
    sigma = std(Rsh);
    p_val = 1 - normcdf(R,mu,sigma);
    reject_h0 = p_val < 0.05/6;
    gam_grid_hd.p_val_shuffle=p_val;
    gam_grid_hd.reject_h0=reject_h0;
    gam_grid_hd.ED_normalized2shuffle=R-mu;
    gam_grid_hd.shuffle_ED_mu=mu;
    gam_grid_hd.shuffle_ED_sigma=sigma;

models_reject=[gam_grid_gz.reject_h0 gam_grid_pl.reject_h0 gam_grid_hd.reject_h0];
models_ED=[gam_grid_gz.ED_normalized2shuffle gam_grid_pl.ED_normalized2shuffle gam_grid_hd.ED_normalized2shuffle];
[m top_model]=max(models_ED)

%%
if sum(models_reject)>0
    
    if sum(models_reject)>1

       switch top_model

           case 1

             [gam_grid_gz_pl, shuffle_gz_pl]=pyGAM_pyr_nogrid_shuffle(y_array,X_arraygz_place,4:6,np,lams,shift);
             
             [gam_grid_gz_hd, shuffle_gz_hd]=pyGAM_pyr_nogrid_shuffle(y_array,X_array_gz_hd,4:6,np,lams,shift);

            temporal.gz_pl{j}=gam_grid_gz_pl;
            temporal.gz_pl_shuffle{j}=shuffle_gz_pl;
            temporal.gz_hd{j}=gam_grid_gz_hd;
            temporal.gz_hd_shuffle{j}=shuffle_gz_hd;
            
            models_reject2=[gam_grid_gz_pl.reject_h0 gam_grid_gz_hd.reject_h0];
            models_ED2=[gam_grid_gz_pl.ED_normalized2shuffle gam_grid_gz_hd.ED_normalized2shuffle];
            [m top_model2]=max(models_ED2)

            if sum(models_reject2)>0

                if sum(models_reject2)>1

                switch top_model2

                    case 1
                      
                        [gam_grid_gz_pl_hd, shuffle_gz_pl_hd]=pyGAM_pyr_nogrid_shuffle(y_array,X_array_gz_place_hd,7:9,np,lams,shift);
                        temporal.gz_pl_hd{j}=gam_grid_gz_pl_hd;
                        temporal.gz_pl_hd_shuffle{j}=shuffle_gz_pl_hd;

                    case 2
                        [gam_grid_gz_pl_hd, shuffle_gz_pl_hd]=pyGAM_pyr_nogrid_shuffle(y_array,X_array_gz_place_hd,4:6,np,lams,shift);
                        temporal.gz_pl_hd{j}=gam_grid_gz_pl_hd;
                        temporal.gz_pl_hd_shuffle{j}=shuffle_gz_pl_hd; 
                end
                else
                rejected=find(models_reject2);
                switch rejected
                    case 1
                [gam_grid_gz_pl_hd, shuffle_gz_pl_hd]=pyGAM_pyr_nogrid_shuffle(y_array,X_array_gz_place_hd,7:9,np,lams,shift);
                temporal.gz_pl_hd{j}=gam_grid_gz_pl_hd;
                temporal.gz_pl_hd_shuffle{j}=shuffle_gz_pl_hd; 

                    case 2
                [gam_grid_gz_pl_hd, shuffle_gz_pl_hd]=pyGAM_pyr_nogrid_shuffle(y_array,X_array_gz_place_hd,4:6,np,lams,shift);
                temporal.gz_pl_hd{j}=gam_grid_gz_pl_hd;
                temporal.gz_pl_hd_shuffle{j}=shuffle_gz_pl_hd; 
                end
                end
            end



                
           
           case 2
             [gam_grid_gz_pl, shuffle_gz_pl]=pyGAM_pyr_nogrid_shuffle(y_array,X_arraygz_place,1:3,np,lams,shift);
             
             [gam_grid_pl_hd, shuffle_pl_hd]=pyGAM_pyr_nogrid_shuffle(y_array,X_array_place_hd,4:6,np,lams,shift);

            temporal.gz_pl{j}=gam_grid_gz_pl;
            temporal.gz_pl_shuffle{j}=shuffle_gz_pl;
            temporal.pl_hd{j}=gam_grid_pl_hd;
            temporal.pl_hd_shuffle{j}=shuffle_pl_hd;
            
            models_reject2=[gam_grid_gz_pl.reject_h0 gam_grid_pl_hd.reject_h0];
            models_ED2=[gam_grid_gz_pl.ED_normalized2shuffle gam_grid_pl_hd.ED_normalized2shuffle];
            [m top_model2]=max(models_ED2)

            if sum(models_reject2)>0

                if sum(models_reject2)>1

                switch top_model2

                    case 1
                      
                        [gam_grid_gz_pl_hd, shuffle_gz_pl_hd]=pyGAM_pyr_nogrid_shuffle(y_array,X_array_gz_place_hd,7:9,np,lams,shift);
                        temporal.gz_pl_hd{j}=gam_grid_gz_pl_hd;
                        temporal.gz_pl_hd_shuffle{j}=shuffle_gz_pl_hd;

                    case 2
                        [gam_grid_gz_pl_hd, shuffle_gz_pl_hd]=pyGAM_pyr_nogrid_shuffle(y_array,X_array_gz_place_hd,1:3,np,lams,shift);
                        temporal.gz_pl_hd{j}=gam_grid_gz_pl_hd;
                        temporal.gz_pl_hd_shuffle{j}=shuffle_gz_pl_hd; 
                end
                else
                rejected=find(models_reject2);
                switch rejected
                    case 1
                [gam_grid_gz_pl_hd, shuffle_gz_pl_hd]=pyGAM_pyr_nogrid_shuffle(y_array,X_array_gz_place_hd,7:9,np,lams,shift);
                temporal.gz_pl_hd{j}=gam_grid_gz_pl_hd;
                temporal.gz_pl_hd_shuffle{j}=shuffle_gz_pl_hd; 

                    case 2
                [gam_grid_gz_pl_hd, shuffle_gz_pl_hd]=pyGAM_pyr_nogrid_shuffle(y_array,X_array_gz_place_hd,1:3,np,lams,shift);
                temporal.gz_pl_hd{j}=gam_grid_gz_pl_hd;
                temporal.gz_pl_hd_shuffle{j}=shuffle_gz_pl_hd; 
                end
                end
            end

          

       case 3
             [gam_grid_gz_hd, shuffle_gz_hd]=pyGAM_pyr_nogrid_shuffle(y_array,X_array_gz_hd,1:3,np,lams,shift);
             
             [gam_grid_pl_hd, shuffle_pl_hd]=pyGAM_pyr_nogrid_shuffle(y_array,X_array_place_hd,1:3,np,lams,shift);

            temporal.gz_hd{j}=gam_grid_gz_hd;
            temporal.gz_hd_shuffle{j}=shuffle_gz_hd;
            temporal.pl_hd{j}=gam_grid_pl_hd;
            temporal.pl_hd_shuffle{j}=shuffle_pl_hd;
            
            models_reject2=[gam_grid_gz_hd.reject_h0 gam_grid_pl_hd.reject_h0];
            models_ED2=[gam_grid_gz_hd.ED_normalized2shuffle gam_grid_pl_hd.ED_normalized2shuffle];
            [m top_model2]=max(models_ED2)

            if sum(models_reject2)>0

                if sum(models_reject2)>1

                switch top_model2

                    case 1
                      
                        [gam_grid_gz_pl_hd, shuffle_gz_pl_hd]=pyGAM_pyr_nogrid_shuffle(y_array,X_array_gz_place_hd,4:6,np,lams,shift);
                        temporal.gz_pl_hd{j}=gam_grid_gz_pl_hd;
                        temporal.gz_pl_hd_shuffle{j}=shuffle_gz_pl_hd;

                    case 2
                        [gam_grid_gz_pl_hd, shuffle_gz_pl_hd]=pyGAM_pyr_nogrid_shuffle(y_array,X_array_gz_place_hd,1:3,np,lams,shift);
                        temporal.gz_pl_hd{j}=gam_grid_gz_pl_hd;
                        temporal.gz_pl_hd_shuffle{j}=shuffle_gz_pl_hd; 
                end
                else
                rejected=find(models_reject2);
                switch rejected
                    case 1
                [gam_grid_gz_pl_hd, shuffle_gz_pl_hd]=pyGAM_pyr_nogrid_shuffle(y_array,X_array_gz_place_hd,4:6,np,lams,shift);
                temporal.gz_pl_hd{j}=gam_grid_gz_pl_hd;
                temporal.gz_pl_hd_shuffle{j}=shuffle_gz_pl_hd; 

                    case 2
                [gam_grid_gz_pl_hd, shuffle_gz_pl_hd]=pyGAM_pyr_nogrid_shuffle(y_array,X_array_gz_place_hd,1:3,np,lams,shift);
                temporal.gz_pl_hd{j}=gam_grid_gz_pl_hd;
                temporal.gz_pl_hd_shuffle{j}=shuffle_gz_pl_hd; 
                end
                end
            end
       end

    else
    rejected3=find(models_reject);
        
    switch rejected3

           case 1

             [gam_grid_gz_pl, shuffle_gz_pl]=pyGAM_pyr_nogrid_shuffle(y_array,X_arraygz_place,4:6,np,lams,shift);
             
             [gam_grid_gz_hd, shuffle_gz_hd]=pyGAM_pyr_nogrid_shuffle(y_array,X_array_gz_hd,4:6,np,lams,shift);

            temporal.gz_pl{j}=gam_grid_gz_pl;
            temporal.gz_pl_shuffle{j}=shuffle_gz_pl;
            temporal.gz_hd{j}=gam_grid_gz_hd;
            temporal.gz_hd_shuffle{j}=shuffle_gz_hd;
            
            models_reject2=[gam_grid_gz_pl.reject_h0 gam_grid_gz_hd.reject_h0];
            models_ED2=[gam_grid_gz_pl.ED_normalized2shuffle gam_grid_gz_hd.ED_normalized2shuffle];
            [m top_model2]=max(models_ED2)

            if sum(models_reject2)>0

                if sum(models_reject2)>1

                switch top_model2

                    case 1
                      
                        [gam_grid_gz_pl_hd, shuffle_gz_pl_hd]=pyGAM_pyr_nogrid_shuffle(y_array,X_array_gz_place_hd,7:9,np,lams,shift);
                        temporal.gz_pl_hd{j}=gam_grid_gz_pl_hd;
                        temporal.gz_pl_hd_shuffle{j}=shuffle_gz_pl_hd;

                    case 2
                        [gam_grid_gz_pl_hd, shuffle_gz_pl_hd]=pyGAM_pyr_nogrid_shuffle(y_array,X_array_gz_place_hd,4:6,np,lams,shift);
                        temporal.gz_pl_hd{j}=gam_grid_gz_pl_hd;
                        temporal.gz_pl_hd_shuffle{j}=shuffle_gz_pl_hd; 
                end
                else
                rejected=find(models_reject2);
                switch rejected
                    case 1
                [gam_grid_gz_pl_hd, shuffle_gz_pl_hd]=pyGAM_pyr_nogrid_shuffle(y_array,X_array_gz_place_hd,7:9,np,lams,shift);
                temporal.gz_pl_hd{j}=gam_grid_gz_pl_hd;
                temporal.gz_pl_hd_shuffle{j}=shuffle_gz_pl_hd; 

                    case 2
                [gam_grid_gz_pl_hd, shuffle_gz_pl_hd]=pyGAM_pyr_nogrid_shuffle(y_array,X_array_gz_place_hd,4:6,np,lams,shift);
                temporal.gz_pl_hd{j}=gam_grid_gz_pl_hd;
                temporal.gz_pl_hd_shuffle{j}=shuffle_gz_pl_hd; 
                end
                end
            end



                
           
           case 2
             [gam_grid_gz_pl, shuffle_gz_pl]=pyGAM_pyr_nogrid_shuffle(y_array,X_arraygz_place,1:3,np,lams,shift);
             
             [gam_grid_pl_hd, shuffle_pl_hd]=pyGAM_pyr_nogrid_shuffle(y_array,X_array_place_hd,4:6,np,lams,shift);

            temporal.gz_pl{j}=gam_grid_gz_pl;
            temporal.gz_pl_shuffle{j}=shuffle_gz_pl;
            temporal.gz_hd{j}=gam_grid_pl_hd;
            temporal.gz_hd_shuffle{j}=shuffle_pl_hd;
            
            models_reject2=[gam_grid_gz_pl.reject_h0 gam_grid_pl_hd.reject_h0];
            models_ED2=[gam_grid_gz_pl.ED_normalized2shuffle gam_grid_pl_hd.ED_normalized2shuffle];
            [m top_model2]=max(models_ED2)

            if sum(models_reject2)>0

                if sum(models_reject2)>1

                switch top_model2

                    case 1
                      
                        [gam_grid_gz_pl_hd, shuffle_gz_pl_hd]=pyGAM_pyr_nogrid_shuffle(y_array,X_array_gz_place_hd,7:9,np,lams,shift);
                        temporal.gz_pl_hd{j}=gam_grid_gz_pl_hd;
                        temporal.gz_pl_hd_shuffle{j}=shuffle_gz_pl_hd;

                    case 2
                        [gam_grid_gz_pl_hd, shuffle_gz_pl_hd]=pyGAM_pyr_nogrid_shuffle(y_array,X_array_gz_place_hd,1:3,np,lams,shift);
                        temporal.gz_pl_hd{j}=gam_grid_gz_pl_hd;
                        temporal.gz_pl_hd_shuffle{j}=shuffle_gz_pl_hd; 
                end
                else
                rejected=find(models_reject2);
                switch rejected
                    case 1
                [gam_grid_gz_pl_hd, shuffle_gz_pl_hd]=pyGAM_pyr_nogrid_shuffle(y_array,X_array_gz_place_hd,7:9,np,lams,shift);
                temporal.gz_pl_hd{j}=gam_grid_gz_pl_hd;
                temporal.gz_pl_hd_shuffle{j}=shuffle_gz_pl_hd; 

                    case 2
                [gam_grid_gz_pl_hd, shuffle_gz_pl_hd]=pyGAM_pyr_nogrid_shuffle(y_array,X_array_gz_place_hd,1:3,np,lams,shift);
                temporal.gz_pl_hd{j}=gam_grid_gz_pl_hd;
                temporal.gz_pl_hd_shuffle{j}=shuffle_gz_pl_hd; 
                end
                end
            end

          

       case 3
             [gam_grid_gz_hd, shuffle_gz_hd]=pyGAM_pyr_nogrid_shuffle(y_array,X_array_gz_hd,1:3,np,lams,shift);
             
             [gam_grid_pl_hd, shuffle_pl_hd]=pyGAM_pyr_nogrid_shuffle(y_array,X_array_place_hd,1:3,np,lams,shift);

            temporal.gz_hd{j}=gam_grid_gz_hd;
            temporal.gz_hd_shuffle{j}=shuffle_gz_hd;
            temporal.pl_hd{j}=gam_grid_pl_hd;
            temporal.pl_hd_shuffle{j}=shuffle_pl_hd;
            
            models_reject2=[gam_grid_gz_hd.reject_h0 gam_grid_pl_hd.reject_h0];
            models_ED2=[gam_grid_gz_hd.ED_normalized2shuffle gam_grid_pl_hd.ED_normalized2shuffle];
            [m top_model2]=max(models_ED2)

            if sum(models_reject2)>0

                if sum(models_reject2)>1

                switch top_model2

                    case 1
                      
                        [gam_grid_gz_pl_hd, shuffle_gz_pl_hd]=pyGAM_pyr_nogrid_shuffle(y_array,X_array_gz_place_hd,4:6,np,lams,shift);
                        temporal.gz_pl_hd{j}=gam_grid_gz_pl_hd;
                        temporal.gz_pl_hd_shuffle{j}=shuffle_gz_pl_hd;

                    case 2
                        [gam_grid_gz_pl_hd, shuffle_gz_pl_hd]=pyGAM_pyr_nogrid_shuffle(y_array,X_array_gz_place_hd,1:3,np,lams,shift);
                        temporal.gz_pl_hd{j}=gam_grid_gz_pl_hd;
                        temporal.gz_pl_hd_shuffle{j}=shuffle_gz_pl_hd; 
                end
                else
                rejected=find(models_reject2);
                switch rejected
                    case 1
                [gam_grid_gz_pl_hd, shuffle_gz_pl_hd]=pyGAM_pyr_nogrid_shuffle(y_array,X_array_gz_place_hd,4:6,np,lams,shift);
                temporal.gz_pl_hd{j}=gam_grid_gz_pl_hd;
                temporal.gz_pl_hd_shuffle{j}=shuffle_gz_pl_hd; 

                    case 2
                [gam_grid_gz_pl_hd, shuffle_gz_pl_hd]=pyGAM_pyr_nogrid_shuffle(y_array,X_array_gz_place_hd,1:3,np,lams,shift);
                temporal.gz_pl_hd{j}=gam_grid_gz_pl_hd;
                temporal.gz_pl_hd_shuffle{j}=shuffle_gz_pl_hd; 
                end
                end
            end
    end
    end

end
%%
temporal.hd_shuffle{j}=gam_grid_hd_shuffle;
temporal.gz_shuffle{j}=gam_grid_gz_shuffle;
temporal.pl_shuffle{j}=gam_grid_pl_shuffle;
temporal.hd{j}=gam_grid_hd;
temporal.gz{j}=gam_grid_gz;
temporal.pl{j}=gam_grid_pl;
else
    continue
end

end

if any(strcmp('gz',fieldnames(temporal)))
pyGAM=[pyGAM;temporal];
end
