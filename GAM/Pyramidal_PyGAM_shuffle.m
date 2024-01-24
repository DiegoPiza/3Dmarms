%% This script runs the GAM analysis described in methods, and replicates figure S5c
% To install GAM_models, follow the instructions here (https://GAM_models.readthedocs.io/en/latest/) for pip installation
% gz models correspond to view, pl to place and hd to head direction
clear all
close all
rng shuffle
load sample_data.mat
variable_types = repmat(cellstr('cell'), 1, 15);
shuffle_number=100;%how many shuffles
grid_number=100;%how many grid search lambda values
neuron='S20190619C27_18' ; %ID of Neuron to model
n_folds=5; % how many crossvalidation partitions
sr=60;% data sampling rate



%Replace bellow with your python environment path
setenv('path',['...\AppData\Local\Programs\Python', getenv('path')]);
pe=pyenv('Version','...\AppData\Local\Programs\Python\Python310\pythonw.exe',...
    'ExecutionMode','InProcess');


%%
types_ID=string(compactWave.UnitIDs);
lams = rand([grid_number 3]); %labda values (smoothing penalization)
lams = lams * 6 - 3 ;
lams = py.numpy.array(10.^lams) ;
TrackingData=even_TData(TrackingData);
puls=TrackingData.neural_time;
fields=fieldnames(units.singleunits);
heatmrawf=Head_Orientation_No_floors(TrackingData); %obtains view
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
GAM_models=table('Size',[1,15],'VariableTypes',variable_types,'VariableNames',variable_names);



unit=units.singleunits.(neuron);

unitID=string(neuron);
type_i=find(unitID==types_ID);
if isempty(type_i)
    error('no neuronal type info')
end

neuron_type=compactWave.NeuronType(type_i);

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

close all
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

X_arrayhd=py.numpy.array([x_r y_r z_r]);
X_arraygz=py.numpy.array([Xgf Ygf Zgf]);
X_arrayplace=py.numpy.array([X Y Z]);

X_arraygz_place=py.numpy.array([Xgf Ygf Zgf X Y Z]);

X_array_gz_place_hd=py.numpy.array([Xgf Ygf Zgf X Y Z x_r y_r z_r]);

X_array_gz_hd=py.numpy.array([Xgf Ygf Zgf x_r y_r z_r]);

X_array_place_hd=py.numpy.array([X Y Z x_r y_r z_r]);

X_array_gz_hd=py.numpy.array([Xgf Ygf Zgf x_r y_r z_r]);


y_array=py.numpy.array(spikes);


if neuron_type=="Pyr"

    if isempty(GAM_models.unitID{1})
        GAM_models.unitID{1}=unitID;
        jj=1;
    else
        jj=height(GAM_models)+1;
        GAM_models.unitID{jj}=unitID;
    end


    gamgz=np.GAM(pyargs('distribution','poisson','link','log'));
    gam_grid_gz=gamgz.gridsearch(X_arraygz ,y_array,pyargs('lam',lams));
    gam_grid_gz=pyGAM2struct(gam_grid_gz);
    lamgz=gam_grid_gz.lam;

    shift=randi([5*sr length(spikes)-(5*sr)],[shuffle_number,1]); %-+ 5 sec shift

    shuffle_gz=cell(shuffle_number,1);
    shuffle_ev_gz=nan(shuffle_number,1);
    for i=1:shuffle_number
        spikes_sh=circshift(spikes,shift(i));
        y_array_gz_sh=py.numpy.array(spikes_sh);
        gam_grid_gz_shuffle=np.GAM(pyargs('distribution','poisson','link','log','lam',lamgz));
        gam_grid_gz_shuffle=gam_grid_gz_shuffle.fit(X_arraygz ,y_array_gz_sh);
        gam_grid_gz_shuffle=pyGAM2struct(gam_grid_gz_shuffle);
        shuffle_gz{i}=gam_grid_gz_shuffle;
        shuffle_ev_gz(i)=gam_grid_gz_shuffle.pseudo_r2.explained_deviance;
    end


    %statistics
    Rsh=shuffle_ev_gz;
    R=gam_grid_gz.pseudo_r2.explained_deviance;
    mu = mean(Rsh);
    if isnan(mu)
        error('nan mean detected')
    end
    sigma = std(Rsh);
    p_val = 1 - normcdf(R,mu,sigma);
    reject_h0 = p_val < 0.05/6;
    gam_grid_gz.p_val_shuffle=p_val;
    gam_grid_gz.reject_h0=reject_h0;
    gam_grid_gz.ED_normalized2shuffle=R-mu;
    gam_grid_gz.shuffle_ED_mu=mu;
    gam_grid_gz.shuffle_ED_sigma=sigma;



    gampl=np.GAM(pyargs('distribution','poisson','link','log'));
    gam_grid_pl=gampl.gridsearch(X_arrayplace ,y_array,pyargs('lam',lams));
    gam_grid_pl=pyGAM2struct(gam_grid_pl);
    lampl=gam_grid_pl.lam;

    shuffle_pl=cell(shuffle_number,1);
    shuffle_ev_pl=nan(shuffle_number,1);
    for i=1:shuffle_number
        spikes_sh=circshift(spikes,shift(i));
        y_array_pl_sh=py.numpy.array(spikes_sh);
        gam_grid_pl_shuffle=np.GAM(pyargs('distribution','poisson','link','log','lam',lampl));
        gam_grid_pl_shuffle=gam_grid_pl_shuffle.fit(X_arrayplace ,y_array_pl_sh);
        gam_grid_pl_shuffle=pyGAM2struct(gam_grid_pl_shuffle);
        shuffle_pl{i}=gam_grid_pl_shuffle;
        shuffle_ev_pl(i)=gam_grid_pl_shuffle.pseudo_r2.explained_deviance;
    end

    %statistics
    Rsh=shuffle_ev_pl;
    R=gam_grid_pl.pseudo_r2.explained_deviance;
    mu = mean(Rsh);
    if isnan(mu)
        error('nan mean detected')
    end
    sigma = std(Rsh);
    p_val = 1 - normcdf(R,mu,sigma);
    reject_h0 = p_val < 0.05/6;
    gam_grid_pl.p_val_shuffle=p_val;
    gam_grid_pl.reject_h0=reject_h0;
    gam_grid_pl.ED_normalized2shuffle=R-mu;
    gam_grid_pl.shuffle_ED_mu=mu;
    gam_grid_pl.shuffle_ED_sigma=sigma;


    gamhd=np.GAM(pyargs('distribution','poisson','link','log'));
    gam_grid_hd=gamhd.gridsearch(X_arrayhd ,y_array,pyargs('lam',lams));
    gam_grid_hd=pyGAM2struct(gam_grid_hd);
    lamhd=gam_grid_hd.lam;

    shuffle_hd=cell(shuffle_number,1);
    shuffle_ev_hd=nan(shuffle_number,1);
    for i=1:shuffle_number
        spikes_sh=circshift(spikes,shift(i));
        y_array_hd_sh=py.numpy.array(spikes_sh);
        gam_grid_hd_shuffle=np.GAM(pyargs('distribution','poisson','link','log','lam',lamhd));
        gam_grid_hd_shuffle=gam_grid_hd_shuffle.fit(X_arrayhd ,y_array_hd_sh);
        gam_grid_hd_shuffle=pyGAM2struct(gam_grid_hd_shuffle);
        shuffle_hd{i}=gam_grid_hd_shuffle;
        shuffle_ev_hd(i)=gam_grid_hd_shuffle.pseudo_r2.explained_deviance;
    end

    %statistics
    Rsh=shuffle_ev_hd;
    R=gam_grid_hd.pseudo_r2.explained_deviance;
    mu = mean(Rsh);
    if isnan(mu)
        error('nan mean detected')
    end
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

                    [gam_grid_gz_pl, shuffle_gz_pl]=pyGAM_int_2models_shuffle(y_array,X_arraygz_place,4:6,np,[lamgz lampl],shift,shuffle_number);

                    [gam_grid_gz_hd, shuffle_gz_hd]=pyGAM_int_2models_shuffle(y_array,X_array_gz_hd,4:6,np,[lamgz lamhd],shift,shuffle_number);

                    GAM_models.gz_pl{jj}=gam_grid_gz_pl;
                    GAM_models.gz_pl_shuffle{jj}=shuffle_gz_pl;
                    GAM_models.gz_hd{jj}=gam_grid_gz_hd;
                    GAM_models.gz_hd_shuffle{jj}=shuffle_gz_hd;

                    models_reject2=[gam_grid_gz_pl.reject_h0 gam_grid_gz_hd.reject_h0];
                    models_ED2=[gam_grid_gz_pl.ED_normalized2shuffle gam_grid_gz_hd.ED_normalized2shuffle];
                    [m top_model2]=max(models_ED2)

                    if sum(models_reject2)>0

                        if sum(models_reject2)>1

                            switch top_model2

                                case 1

                                    [gam_grid_gz_pl_hd, shuffle_gz_pl_hd]=pyGAM_int_2models_shuffle(y_array,X_array_gz_place_hd,7:9,np,[lamgz lampl lamhd],shift,shuffle_number);
                                    GAM_models.gz_pl_hd{jj}=gam_grid_gz_pl_hd;
                                    GAM_models.gz_pl_hd_shuffle{jj}=shuffle_gz_pl_hd;

                                case 2
                                    [gam_grid_gz_pl_hd, shuffle_gz_pl_hd]=pyGAM_int_2models_shuffle(y_array,X_array_gz_place_hd,4:6,np,[lamgz lampl lamhd],shift,shuffle_number);
                                    GAM_models.gz_pl_hd{jj}=gam_grid_gz_pl_hd;
                                    GAM_models.gz_pl_hd_shuffle{jj}=shuffle_gz_pl_hd;
                            end
                        else
                            rejected=find(models_reject2);
                            switch rejected
                                case 1
                                    [gam_grid_gz_pl_hd, shuffle_gz_pl_hd]=pyGAM_int_2models_shuffle(y_array,X_array_gz_place_hd,7:9,np,[lamgz lampl lamhd],shift,shuffle_number);
                                    GAM_models.gz_pl_hd{jj}=gam_grid_gz_pl_hd;
                                    GAM_models.gz_pl_hd_shuffle{jj}=shuffle_gz_pl_hd;

                                case 2
                                    [gam_grid_gz_pl_hd, shuffle_gz_pl_hd]=pyGAM_int_2models_shuffle(y_array,X_array_gz_place_hd,4:6,np,[lamgz lampl lamhd],shift,shuffle_number);
                                    GAM_models.gz_pl_hd{jj}=gam_grid_gz_pl_hd;
                                    GAM_models.gz_pl_hd_shuffle{jj}=shuffle_gz_pl_hd;
                            end
                        end
                    end





                case 2
                    [gam_grid_gz_pl, shuffle_gz_pl]=pyGAM_int_2models_shuffle(y_array,X_arraygz_place,1:3,np,[lamgz lampl],shift,shuffle_number);

                    [gam_grid_pl_hd, shuffle_pl_hd]=pyGAM_int_2models_shuffle(y_array,X_array_place_hd,4:6,np,[lampl lamhd],shift,shuffle_number);

                    GAM_models.gz_pl{jj}=gam_grid_gz_pl;
                    GAM_models.gz_pl_shuffle{jj}=shuffle_gz_pl;
                    GAM_models.pl_hd{jj}=gam_grid_pl_hd;
                    GAM_models.pl_hd_shuffle{jj}=shuffle_pl_hd;

                    models_reject2=[gam_grid_gz_pl.reject_h0 gam_grid_pl_hd.reject_h0];
                    models_ED2=[gam_grid_gz_pl.ED_normalized2shuffle gam_grid_pl_hd.ED_normalized2shuffle];
                    [m top_model2]=max(models_ED2)

                    if sum(models_reject2)>0

                        if sum(models_reject2)>1

                            switch top_model2

                                case 1

                                    [gam_grid_gz_pl_hd, shuffle_gz_pl_hd]=pyGAM_int_2models_shuffle(y_array,X_array_gz_place_hd,7:9,np,[lamgz lampl lamhd],shift,shuffle_number);
                                    GAM_models.gz_pl_hd{jj}=gam_grid_gz_pl_hd;
                                    GAM_models.gz_pl_hd_shuffle{jj}=shuffle_gz_pl_hd;

                                case 2
                                    [gam_grid_gz_pl_hd, shuffle_gz_pl_hd]=pyGAM_int_2models_shuffle(y_array,X_array_gz_place_hd,1:3,np,[lamgz lampl lamhd],shift,shuffle_number);
                                    GAM_models.gz_pl_hd{jj}=gam_grid_gz_pl_hd;
                                    GAM_models.gz_pl_hd_shuffle{jj}=shuffle_gz_pl_hd;
                            end
                        else
                            rejected=find(models_reject2);
                            switch rejected
                                case 1
                                    [gam_grid_gz_pl_hd, shuffle_gz_pl_hd]=pyGAM_int_2models_shuffle(y_array,X_array_gz_place_hd,7:9,np,[lamgz lampl lamhd],shift,shuffle_number);
                                    GAM_models.gz_pl_hd{jj}=gam_grid_gz_pl_hd;
                                    GAM_models.gz_pl_hd_shuffle{jj}=shuffle_gz_pl_hd;

                                case 2
                                    [gam_grid_gz_pl_hd, shuffle_gz_pl_hd]=pyGAM_int_2models_shuffle(y_array,X_array_gz_place_hd,1:3,np,[lamgz lampl lamhd],shift,shuffle_number);
                                    GAM_models.gz_pl_hd{jj}=gam_grid_gz_pl_hd;
                                    GAM_models.gz_pl_hd_shuffle{jj}=shuffle_gz_pl_hd;
                            end
                        end
                    end



                case 3
                    [gam_grid_gz_hd, shuffle_gz_hd]=pyGAM_int_2models_shuffle(y_array,X_array_gz_hd,1:3,np,[lamgz lamhd],shift,shuffle_number);

                    [gam_grid_pl_hd, shuffle_pl_hd]=pyGAM_int_2models_shuffle(y_array,X_array_place_hd,1:3,np,[lampl lamhd],shift,shuffle_number);

                    GAM_models.gz_hd{jj}=gam_grid_gz_hd;
                    GAM_models.gz_hd_shuffle{jj}=shuffle_gz_hd;
                    GAM_models.pl_hd{jj}=gam_grid_pl_hd;
                    GAM_models.pl_hd_shuffle{jj}=shuffle_pl_hd;

                    models_reject2=[gam_grid_gz_hd.reject_h0 gam_grid_pl_hd.reject_h0];
                    models_ED2=[gam_grid_gz_hd.ED_normalized2shuffle gam_grid_pl_hd.ED_normalized2shuffle];
                    [m top_model2]=max(models_ED2)

                    if sum(models_reject2)>0

                        if sum(models_reject2)>1

                            switch top_model2

                                case 1

                                    [gam_grid_gz_pl_hd, shuffle_gz_pl_hd]=pyGAM_int_2models_shuffle(y_array,X_array_gz_place_hd,4:6,np,[lamgz lampl lamhd],shift,shuffle_number);
                                    GAM_models.gz_pl_hd{jj}=gam_grid_gz_pl_hd;
                                    GAM_models.gz_pl_hd_shuffle{jj}=shuffle_gz_pl_hd;

                                case 2
                                    [gam_grid_gz_pl_hd, shuffle_gz_pl_hd]=pyGAM_int_2models_shuffle(y_array,X_array_gz_place_hd,1:3,np,[lamgz lampl lamhd],shift,shuffle_number);
                                    GAM_models.gz_pl_hd{jj}=gam_grid_gz_pl_hd;
                                    GAM_models.gz_pl_hd_shuffle{jj}=shuffle_gz_pl_hd;
                            end
                        else
                            rejected=find(models_reject2);
                            switch rejected
                                case 1
                                    [gam_grid_gz_pl_hd, shuffle_gz_pl_hd]=pyGAM_int_2models_shuffle(y_array,X_array_gz_place_hd,4:6,np,[lamgz lampl lamhd],shift,shuffle_number);
                                    GAM_models.gz_pl_hd{jj}=gam_grid_gz_pl_hd;
                                    GAM_models.gz_pl_hd_shuffle{jj}=shuffle_gz_pl_hd;

                                case 2
                                    [gam_grid_gz_pl_hd, shuffle_gz_pl_hd]=pyGAM_int_2models_shuffle(y_array,X_array_gz_place_hd,1:3,np,[lamgz lampl lamhd],shift,shuffle_number);
                                    GAM_models.gz_pl_hd{jj}=gam_grid_gz_pl_hd;
                                    GAM_models.gz_pl_hd_shuffle{jj}=shuffle_gz_pl_hd;
                            end
                        end
                    end
            end

        else
            rejected3=find(models_reject);

            switch rejected3

                case 1

                    [gam_grid_gz_pl, shuffle_gz_pl]=pyGAM_int_2models_shuffle(y_array,X_arraygz_place,4:6,np,[lamgz lampl],shift,shuffle_number);

                    [gam_grid_gz_hd, shuffle_gz_hd]=pyGAM_int_2models_shuffle(y_array,X_array_gz_hd,4:6,np,[lamgz lamhd],shift,shuffle_number);

                    GAM_models.gz_pl{jj}=gam_grid_gz_pl;
                    GAM_models.gz_pl_shuffle{jj}=shuffle_gz_pl;
                    GAM_models.gz_hd{jj}=gam_grid_gz_hd;
                    GAM_models.gz_hd_shuffle{jj}=shuffle_gz_hd;

                    models_reject2=[gam_grid_gz_pl.reject_h0 gam_grid_gz_hd.reject_h0];
                    models_ED2=[gam_grid_gz_pl.ED_normalized2shuffle gam_grid_gz_hd.ED_normalized2shuffle];
                    [m top_model2]=max(models_ED2)

                    if sum(models_reject2)>0

                        if sum(models_reject2)>1

                            switch top_model2

                                case 1

                                    [gam_grid_gz_pl_hd, shuffle_gz_pl_hd]=pyGAM_int_2models_shuffle(y_array,X_array_gz_place_hd,7:9,np,[lamgz lampl lamhd],shift,shuffle_number);
                                    GAM_models.gz_pl_hd{jj}=gam_grid_gz_pl_hd;
                                    GAM_models.gz_pl_hd_shuffle{jj}=shuffle_gz_pl_hd;

                                case 2
                                    [gam_grid_gz_pl_hd, shuffle_gz_pl_hd]=pyGAM_int_2models_shuffle(y_array,X_array_gz_place_hd,4:6,np,[lamgz lampl lamhd],shift,shuffle_number);
                                    GAM_models.gz_pl_hd{jj}=gam_grid_gz_pl_hd;
                                    GAM_models.gz_pl_hd_shuffle{jj}=shuffle_gz_pl_hd;
                            end
                        else
                            rejected=find(models_reject2);
                            switch rejected
                                case 1
                                    [gam_grid_gz_pl_hd, shuffle_gz_pl_hd]=pyGAM_int_2models_shuffle(y_array,X_array_gz_place_hd,7:9,np,[lamgz lampl lamhd],shift,shuffle_number);
                                    GAM_models.gz_pl_hd{jj}=gam_grid_gz_pl_hd;
                                    GAM_models.gz_pl_hd_shuffle{jj}=shuffle_gz_pl_hd;

                                case 2
                                    [gam_grid_gz_pl_hd, shuffle_gz_pl_hd]=pyGAM_int_2models_shuffle(y_array,X_array_gz_place_hd,4:6,np,[lamgz lampl lamhd],shift,shuffle_number);
                                    GAM_models.gz_pl_hd{jj}=gam_grid_gz_pl_hd;
                                    GAM_models.gz_pl_hd_shuffle{jj}=shuffle_gz_pl_hd;
                            end
                        end
                    end





                case 2
                    [gam_grid_gz_pl, shuffle_gz_pl]=pyGAM_int_2models_shuffle(y_array,X_arraygz_place,1:3,np,[lamgz lampl],shift,shuffle_number);

                    [gam_grid_pl_hd, shuffle_pl_hd]=pyGAM_int_2models_shuffle(y_array,X_array_place_hd,4:6,np,[lampl lamhd],shift,shuffle_number);

                    GAM_models.gz_pl{jj}=gam_grid_gz_pl;
                    GAM_models.gz_pl_shuffle{jj}=shuffle_gz_pl;
                    GAM_models.gz_hd{jj}=gam_grid_pl_hd;
                    GAM_models.gz_hd_shuffle{jj}=shuffle_pl_hd;

                    models_reject2=[gam_grid_gz_pl.reject_h0 gam_grid_pl_hd.reject_h0];
                    models_ED2=[gam_grid_gz_pl.ED_normalized2shuffle gam_grid_pl_hd.ED_normalized2shuffle];
                    [m top_model2]=max(models_ED2)

                    if sum(models_reject2)>0

                        if sum(models_reject2)>1

                            switch top_model2

                                case 1

                                    [gam_grid_gz_pl_hd, shuffle_gz_pl_hd]=pyGAM_int_2models_shuffle(y_array,X_array_gz_place_hd,7:9,np,[lamgz lampl lamhd],shift,shuffle_number);
                                    GAM_models.gz_pl_hd{jj}=gam_grid_gz_pl_hd;
                                    GAM_models.gz_pl_hd_shuffle{jj}=shuffle_gz_pl_hd;

                                case 2
                                    [gam_grid_gz_pl_hd, shuffle_gz_pl_hd]=pyGAM_int_2models_shuffle(y_array,X_array_gz_place_hd,1:3,np,[lamgz lampl lamhd],shift,shuffle_number);
                                    GAM_models.gz_pl_hd{jj}=gam_grid_gz_pl_hd;
                                    GAM_models.gz_pl_hd_shuffle{jj}=shuffle_gz_pl_hd;
                            end
                        else
                            rejected=find(models_reject2);
                            switch rejected
                                case 1
                                    [gam_grid_gz_pl_hd, shuffle_gz_pl_hd]=pyGAM_int_2models_shuffle(y_array,X_array_gz_place_hd,7:9,np,[lamgz lampl lamhd],shift,shuffle_number);
                                    GAM_models.gz_pl_hd{jj}=gam_grid_gz_pl_hd;
                                    GAM_models.gz_pl_hd_shuffle{jj}=shuffle_gz_pl_hd;

                                case 2
                                    [gam_grid_gz_pl_hd, shuffle_gz_pl_hd]=pyGAM_int_2models_shuffle(y_array,X_array_gz_place_hd,1:3,np,[lamgz lampl lamhd],shift,shuffle_number);
                                    GAM_models.gz_pl_hd{jj}=gam_grid_gz_pl_hd;
                                    GAM_models.gz_pl_hd_shuffle{jj}=shuffle_gz_pl_hd;
                            end
                        end
                    end



                case 3
                    [gam_grid_gz_hd, shuffle_gz_hd]=pyGAM_int_2models_shuffle(y_array,X_array_gz_hd,1:3,np,[lamgz lamhd],shift,shuffle_number);

                    [gam_grid_pl_hd, shuffle_pl_hd]=pyGAM_int_2models_shuffle(y_array,X_array_place_hd,1:3,np,[lampl lamhd],shift,shuffle_number);

                    GAM_models.gz_hd{jj}=gam_grid_gz_hd;
                    GAM_models.gz_hd_shuffle{jj}=shuffle_gz_hd;
                    GAM_models.pl_hd{jj}=gam_grid_pl_hd;
                    GAM_models.pl_hd_shuffle{jj}=shuffle_pl_hd;

                    models_reject2=[gam_grid_gz_hd.reject_h0 gam_grid_pl_hd.reject_h0];
                    models_ED2=[gam_grid_gz_hd.ED_normalized2shuffle gam_grid_pl_hd.ED_normalized2shuffle];
                    [m top_model2]=max(models_ED2)

                    if sum(models_reject2)>0

                        if sum(models_reject2)>1

                            switch top_model2

                                case 1

                                    [gam_grid_gz_pl_hd, shuffle_gz_pl_hd]=pyGAM_int_2models_shuffle(y_array,X_array_gz_place_hd,4:6,np,[lamgz lampl lamhd],shift,shuffle_number);
                                    GAM_models.gz_pl_hd{jj}=gam_grid_gz_pl_hd;
                                    GAM_models.gz_pl_hd_shuffle{jj}=shuffle_gz_pl_hd;

                                case 2
                                    [gam_grid_gz_pl_hd, shuffle_gz_pl_hd]=pyGAM_int_2models_shuffle(y_array,X_array_gz_place_hd,1:3,np,[lamgz lampl lamhd],shift,shuffle_number);
                                    GAM_models.gz_pl_hd{jj}=gam_grid_gz_pl_hd;
                                    GAM_models.gz_pl_hd_shuffle{jj}=shuffle_gz_pl_hd;
                            end
                        else
                            rejected=find(models_reject2);
                            switch rejected
                                case 1
                                    [gam_grid_gz_pl_hd, shuffle_gz_pl_hd]=pyGAM_int_2models_shuffle(y_array,X_array_gz_place_hd,4:6,np,[lamgz lampl lamhd],shift,shuffle_number);
                                    GAM_models.gz_pl_hd{jj}=gam_grid_gz_pl_hd;
                                    GAM_models.gz_pl_hd_shuffle{jj}=shuffle_gz_pl_hd;

                                case 2
                                    [gam_grid_gz_pl_hd, shuffle_gz_pl_hd]=pyGAM_int_2models_shuffle(y_array,X_array_gz_place_hd,1:3,np,[lamgz lampl lamhd],shift,shuffle_number);
                                    GAM_models.gz_pl_hd{jj}=gam_grid_gz_pl_hd;
                                    GAM_models.gz_pl_hd_shuffle{jj}=shuffle_gz_pl_hd;
                            end
                        end
                    end
            end
        end

    end
    %%
    GAM_models.hd_shuffle{jj}=shuffle_hd;
    GAM_models.gz_shuffle{jj}=shuffle_gz;
    GAM_models.pl_shuffle{jj}=shuffle_pl;
    GAM_models.hd{jj}=gam_grid_hd;
    GAM_models.gz{jj}=gam_grid_gz;
    GAM_models.pl{jj}=gam_grid_pl;

else
    error('neuron is not a pyramidal cell')
end

%% Checking the best model
% bestmodel combinations 1=gz 2=pl 3=hd 4=gzpl 5=gzhd 6= plhd 7= gzplhd

i=height(GAM_models);
structGAM=GAM_models;
hd_dev_exp=[structGAM.hd{i,1}.ED_normalized2shuffle];
gz_dev_exp=[structGAM.gz{i,1}.ED_normalized2shuffle];
pl_dev_exp=[structGAM.pl{i,1}.ED_normalized2shuffle];

hd_rej=[structGAM.hd{i,1}.reject_h0];
gz_rej=[structGAM.gz{i,1}.reject_h0];
pl_rej=[structGAM.pl{i,1}.reject_h0];

if ~isempty(structGAM.gz_pl_hd{i,1})
    secondmodelreject=true;
    thirdmodelreject=structGAM.gz_pl_hd{i,1}.reject_h0;

    %what were the best two variables of the 2nd model
    % 1 = gz+pl 2 = gz+hd 3 = pl+hd
    shuffledpred=[structGAM.gz_pl_hd{i,1}.predictor2shuffle];
    switch shuffledpred(1)

        case 1
            secondmodelidentity=3;

        case 4
            secondmodelidentity=2;

        case 7
            secondmodelidentity=1;
    end

else
    secondmodelreject=[secondmodelreject;false];
    secondmodelidentity=[secondmodelidentity;0];
    thirdmodelreject=[thirdmodelreject;false];
end

dev_exp=[gz_dev_exp pl_dev_exp hd_dev_exp];
[m max_model]=max(dev_exp,[],2);
bestmodelidx=[];
firstmodelreject=[gz_rej pl_rej hd_rej];
firstmodelsig=max_model((sum(firstmodelreject,2)>0)'...
    &~(secondmodelreject>0)'&~(thirdmodelreject>0)');
firstmodelidx=find((sum(firstmodelreject,2)>0)'...
    &~(secondmodelreject>0)'&~(thirdmodelreject>0)');
bestmodelidx(firstmodelidx)=firstmodelsig;
secondmodelsig=secondmodelidentity((sum(firstmodelreject,2)>0)'&(secondmodelreject>0)'&~(thirdmodelreject>0)');
secondmodelidx=find((sum(firstmodelreject,2)>0)'&(secondmodelreject>0)'&~(thirdmodelreject>0)');
% 1 = gz+pl 2 = gz+hd 3 = pl+hd
bestmodelidx(secondmodelidx(secondmodelsig==1))=4;
bestmodelidx(secondmodelidx(secondmodelsig==2))=5;
bestmodelidx(secondmodelidx(secondmodelsig==3))=6;
thirdmodelsig=max_model((sum(firstmodelreject,2)>0)'&(thirdmodelreject>0)');
thirdmodelidx=find((sum(firstmodelreject,2)>0)'&(thirdmodelreject>0)');
bestmodelidx(thirdmodelidx)=7;



%% Stratified crossvalidation of the best model

bestmodel=bestmodelidx;
gamid=height(GAM_models);
j=height(GAM_models);
switch bestmodel
    case 1
        lams=GAM_models.gz{gamid,:}.lam;
        xarray=[Xgf Ygf Zgf];
        [test_indices train_indices ind]=cv_stratified_multilabel_DBP(xarray,n_folds,bestmodel); %for stratified cv

        spikes_model_total={};
        spikes_real_total={};
        modelcorr=nan(1,n_folds);


        for nfold=1:n_folds
            y_array_agv=py.numpy.array(spikes(train_indices(:,nfold)));
            X_arrayang_vel=py.numpy.array([xarray(train_indices(:,nfold),:)]);
            X_arrayang_vel_test=py.numpy.array([xarray(test_indices(:,nfold),:)]);
            spikes_test=spikes(test_indices(:,nfold));
            spikes_real_total{nfold}=spikes_test;
            gamagv=np.GAM(pyargs('distribution','poisson','link','log','lam',lams));
            gamagv=gamagv.fit(X_arrayang_vel ,y_array_agv);
            spikes_model=double(gamagv.predict(X_arrayang_vel_test));
            xr=corr(spikes_model',spikes_test);
            spikes_model_total{nfold}=spikes_model;
            modelcorr(nfold)=xr;
        end
    case 2
        lams=GAM_models.pl{gamid,:}.lam;
        xarray=[X Y Z];

        [test_indices train_indices ind]=cv_stratified_multilabel_DBP(xarray,n_folds,bestmodel); %for stratified cv
        spikes_model_total={};
        spikes_real_total={};
        modelcorr=nan(1,n_folds);
        for nfold=1:n_folds
            y_array_agv=py.numpy.array(spikes(train_indices(:,nfold)));
            X_arrayang_vel=py.numpy.array([xarray(train_indices(:,nfold),:)]);
            X_arrayang_vel_test=py.numpy.array([xarray(test_indices(:,nfold),:)]);
            spikes_test=spikes(test_indices(:,nfold));
            spikes_real_total{nfold}=spikes_test;
            gamagv=np.GAM(pyargs('distribution','poisson','link','log','lam',lams));
            gamagv=gamagv.fit(X_arrayang_vel ,y_array_agv);
            spikes_model=double(gamagv.predict(X_arrayang_vel_test));
            xr=corr(spikes_model',spikes_test);
            spikes_model_total{nfold}=spikes_model;
            modelcorr(nfold)=xr;
        end
    case 3
        lams=GAM_models.hd{gamid,:}.lam;
        xarray=[x_r y_r z_r];

        [test_indices train_indices ind]=cv_stratified_multilabel_DBP(xarray,n_folds,bestmodel); %for stratified cv
        spikes_model_total={};
        spikes_real_total={};
        modelcorr=nan(1,n_folds);
        for nfold=1:n_folds
            y_array_agv=py.numpy.array(spikes(train_indices(:,nfold)));
            X_arrayang_vel=py.numpy.array([xarray(train_indices(:,nfold),:)]);
            X_arrayang_vel_test=py.numpy.array([xarray(test_indices(:,nfold),:)]);
            spikes_test=spikes(test_indices(:,nfold));
            spikes_real_total{nfold}=spikes_test;
            gamagv=np.GAM(pyargs('distribution','poisson','link','log','lam',lams));
            gamagv=gamagv.fit(X_arrayang_vel ,y_array_agv);
            spikes_model=double(gamagv.predict(X_arrayang_vel_test));
            xr=corr(spikes_model',spikes_test);
            spikes_model_total{nfold}=spikes_model;
            modelcorr(nfold)=xr;
        end
    case 4
        lams=GAM_models.gz_pl{gamid,:}.lam;
        xarray=[Xgf Ygf Zgf X Y Z];

        [test_indices train_indices ind]=cv_stratified_multilabel_DBP(xarray,n_folds,bestmodel); %for stratified cv
        spikes_model_total={};
        spikes_real_total={};
        modelcorr=nan(1,n_folds);
        for nfold=1:n_folds
            y_array_agv=py.numpy.array(spikes(train_indices(:,nfold)));
            X_arrayang_vel=py.numpy.array([xarray(train_indices(:,nfold),:)]);
            X_arrayang_vel_test=py.numpy.array([xarray(test_indices(:,nfold),:)]);
            spikes_test=spikes(test_indices(:,nfold));
            spikes_real_total{nfold}=spikes_test;
            gamagv=np.GAM(pyargs('distribution','poisson','link','log','lam',lams));
            gamagv=gamagv.fit(X_arrayang_vel ,y_array_agv);
            spikes_model=double(gamagv.predict(X_arrayang_vel_test));
            xr=corr(spikes_model',spikes_test);
            spikes_model_total{nfold}=spikes_model;
            modelcorr(nfold)=xr;
        end
    case 5
        lams=GAM_models.gz_hd{gamid,:}.lam;
        xarray=[Xgf Ygf Zgf x_r y_r z_r];

        [test_indices train_indices ind]=cv_stratified_multilabel_DBP(xarray,n_folds,bestmodel); %for stratified cv
        spikes_model_total={};
        spikes_real_total={};
        modelcorr=nan(1,n_folds);
        for nfold=1:n_folds
            y_array_agv=py.numpy.array(spikes(train_indices(:,nfold)));
            X_arrayang_vel=py.numpy.array([xarray(train_indices(:,nfold),:)]);
            X_arrayang_vel_test=py.numpy.array([xarray(test_indices(:,nfold),:)]);
            spikes_test=spikes(test_indices(:,nfold));
            spikes_real_total{nfold}=spikes_test;
            gamagv=np.GAM(pyargs('distribution','poisson','link','log','lam',lams));
            gamagv=gamagv.fit(X_arrayang_vel ,y_array_agv);
            spikes_model=double(gamagv.predict(X_arrayang_vel_test));
            xr=corr(spikes_model',spikes_test);
            spikes_model_total{nfold}=spikes_model;
            modelcorr(nfold)=xr;
        end

    case 6
        lams=GAM_models.pl_hd{gamid,:}.lam;
        xarray=[X Y Z x_r y_r z_r];

        [test_indices train_indices ind]=cv_stratified_multilabel_DBP(xarray,n_folds,bestmodel); %for stratified cv
        spikes_model_total={};
        spikes_real_total={};
        modelcorr=nan(1,n_folds);
        for nfold=1:n_folds
            y_array_agv=py.numpy.array(spikes(train_indices(:,nfold)));
            X_arrayang_vel=py.numpy.array([xarray(train_indices(:,nfold),:)]);
            X_arrayang_vel_test=py.numpy.array([xarray(test_indices(:,nfold),:)]);
            spikes_test=spikes(test_indices(:,nfold));
            spikes_real_total{nfold}=spikes_test;
            gamagv=np.GAM(pyargs('distribution','poisson','link','log','lam',lams));
            gamagv=gamagv.fit(X_arrayang_vel ,y_array_agv);
            spikes_model=double(gamagv.predict(X_arrayang_vel_test));
            xr=corr(spikes_model',spikes_test);
            spikes_model_total{nfold}=spikes_model;
            modelcorr(nfold)=xr;
        end
    case 7
        lams=GAM_models.gz_pl_hd{gamid,:}.lam;
        xarray=[Xgf Ygf Zgf X Y Z x_r y_r z_r];

        [test_indices train_indices ind]=cv_stratified_multilabel_DBP(xarray,n_folds,bestmodel); %for stratified cv
        spikes_model_total={};
        spikes_real_total={};
        modelcorr=nan(1,n_folds);
        for nfold=1:n_folds
            y_array_agv=py.numpy.array(spikes(train_indices(:,nfold)));
            X_arrayang_vel=py.numpy.array([xarray(train_indices(:,nfold),:)]);
            X_arrayang_vel_test=py.numpy.array([xarray(test_indices(:,nfold),:)]);
            spikes_test=spikes(test_indices(:,nfold));
            spikes_real_total{nfold}=spikes_test;
            gamagv=np.GAM(pyargs('distribution','poisson','link','log','lam',lams));
            gamagv=gamagv.fit(X_arrayang_vel ,y_array_agv);
            spikes_model=double(gamagv.predict(X_arrayang_vel_test));
            xr=corr(spikes_model',spikes_test);
            spikes_model_total{nfold}=spikes_model;
            modelcorr(nfold)=xr;
        end

end

pyGAMcrossval.modelcorr{j}=modelcorr;
pyGAMcrossval.spikes_model_total{j}=spikes_model_total;
pyGAMcrossval.spikes_real_total{j}=spikes_real_total;
pyGAMcrossval.meancorr{j}=mean(modelcorr);
pyGAMcrossval.ind{j}=ind;
pyGAMcrossval.test_indices{j}=test_indices;
pyGAMcrossval.train_indices{j}=train_indices;

%% R squared F-test

p=[];
rejectpyr=[];
R2total=[];
pyr_selec_categ=string;
model_order=string;
for i=1:height(GAM_models)
    % % R2 F-test
    ID=neuron;
    if bestmodelidx(i)==0
        pyr_selec_categ(i)="none";
        p(i)=nan;
        rejectpyr(i)=0;
        model_order(i)="none";
        continue
    end
    R2=mean(pyGAMcrossval.modelcorr{j,1}).^2;
    switch bestmodelidx(i)

        case 1
            k=floor(GAM_models.gz{j,1}.edof);
            pyr_selec_categ(i)="gz";
            model_order(i)="first";
        case 2
            k=floor(GAM_models.pl{j,1}.edof);
            pyr_selec_categ(i)="pl";
            model_order(i)="first";
        case 3
            k=floor(GAM_models.hd{j,1}.edof);
            pyr_selec_categ(i)="hd";
            model_order(i)="first";
        case 4
            k=floor(GAM_models.gz_pl{j,1}.edof);
            pyr_selec_categ(i)="gz_pl";
            model_order(i)="second";
        case 5
            k=floor(GAM_models.gz_hd{j,1}.edof);
            pyr_selec_categ(i)="gz_hd";
            model_order(i)="second";
        case 6
            k=floor(GAM_models.pl_hd{j,1}.edof);
            pyr_selec_categ(i)="pl_hd";
            model_order(i)="second";
        case 7
            k=floor(GAM_models.gz_pl_hd{j,1}.edof);
            pyr_selec_categ(i)="gz_pl_hd";
            model_order(i)="third";
    end
    n=(length(pyGAMcrossval.spikes_real_total{j,1}{1,1})); %https://forum.bionicturtle.com/threads/f-statistic-formula-variations.21735/#:~:text=That%20is%2C%20the%20overall%20F,%3D%20n%20%2D%20k%20%2D%201.
    F = (R2/k) / ((1-R2)/(n-k-1));
    alpha = 0.05;
    Fcrit = finv(1-alpha,k,n-k-1);%https://math.stackexchange.com/questions/617735/multiple-regression-degrees-of-freedom-f-test#:~:text=The%20correct%20approach%20is%20to,freedom%2C%20i.e.%20n%E2%88%921
    p(i) = 1 - fcdf(F,k,n-k-1);
    if p(i) < alpha
        rejectpyr(i)=1
    else
        rejectpyr(i)=0
    end
end
%% Plotting models predictions against firing rate
figure
plot((1:length(spikes))/sr,movmean(spikes,sr/5).*sr,'Color',[0.4 0.4 0.4],'LineWidth',2);
hold on
spikes_model=double(gamagv.predict(X_array_gz_place_hd));
plot((1:length(spikes))/sr,spikes_model.*sr,'Color',[1 0.2 0.2],'LineWidth',1.5);
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





