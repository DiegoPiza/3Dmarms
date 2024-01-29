function structGAM=pyGAM2struct(pGAM)
stats=struct(pGAM.statistics_);
stats.pseudo_r2=struct(stats.pseudo_r2);
% stats.cov=double(stats.cov);
stats.GCV=[];
stats.cov=[];
stats.edof_per_coef=[];%struct(stats.edof_per_coef);
stats.m_features=int64(stats.m_features);
stats.n_samples=int64(stats.n_samples);
stats.se=double(stats.se);
stats.p_values=double(py.numpy.array(stats.p_values));
stats.coef_=double(pGAM.coef_);
terms=struct(pGAM.terms.info);
terms=cell(terms.terms);
for i=1:length(terms)-1
    termst=struct(terms{1,i});
   stats.basis(i)=  string(termst.basis);
   stats.n_splines(i)=   double(termst.n_splines);
   stats.spline_order(i)=   double(termst.spline_order);
   stats.lam(i)=   double(termst.lam);
   stats.dtype(i)=   string(termst.dtype);
end
structGAM=stats;
% if time_diff>duration(1,0,0) 
%    send_email('script update','script still runs','dbuitra2@gmail.com')
%    time=datetime(now,'ConvertFrom','datenum');
% end
end
%%
% structGAM=table;
% structGAM.unitID=pyGAM.unitID;
% 
% for i=1:height(pyGAM) 
% 
%     if any(strcmp('gz_pl_hd',fieldnames(pyGAM)))
%         stats=pyGAM.gz_pl_hd{i,1};
%         if isempty(stats)
%         continue
%         end
%     stats=struct(stats.statistics_);
%     stats.pseudo_r2=struct(stats.pseudo_r2);
%     stats.cov=double(stats.cov);
%     stats.p_values=struct(stats.p_values);
%     stats.edof_per_coef=struct(stats.edof_per_coef);
%     stats.m_features=int64(stats.m_features);
%     stats.n_samples=int64(stats.n_samples);
%     stats.se=double(stats.se);
%     structGAM.gz_pl_hd{i}=stats;
%     end
%     if any(strcmp('hd',fieldnames(pyGAM)))
%         stats=pyGAM.hd{i,1};
%         if isempty(stats)
%         continue
%         end
%     stats=struct(stats.statistics_);
%     stats.pseudo_r2=struct(stats.pseudo_r2);
%     stats.cov=double(stats.cov);
%     stats.p_values=struct(stats.p_values);
%     stats.edof_per_coef=struct(stats.edof_per_coef);
%     stats.m_features=int64(stats.m_features);
%     stats.n_samples=int64(stats.n_samples);
%     stats.se=double(stats.se);
%     structGAM.hd{i}=stats;
%     end
% 
% 
%     if any(strcmp('m1',fieldnames(pyGAM)))
% 
%     stats=pyGAM.m1{i,1};
%     if isempty(stats)
%         continue
%     end
%     stats=struct(stats.statistics_);
%     stats.pseudo_r2=struct(stats.pseudo_r2);
%     stats.cov=double(stats.cov);
%     stats.p_values=struct(stats.p_values);
%     stats.edof_per_coef=struct(stats.edof_per_coef);
%     stats.m_features=int64(stats.m_features);
%     stats.n_samples=int64(stats.n_samples);
%     stats.se=double(stats.se);
%     structGAM.m1{i}=stats;
% 
%     stats=pyGAM.m2{i,1};
%     stats=struct(stats.statistics_);
%     stats.pseudo_r2=struct(stats.pseudo_r2);
%     stats.cov=double(stats.cov);
%     stats.p_values=struct(stats.p_values);
%     stats.edof_per_coef=struct(stats.edof_per_coef);
%     stats.m_features=int64(stats.m_features);
%     stats.n_samples=int64(stats.n_samples);
%     stats.se=double(stats.se);
%     structGAM.m2{i}=stats;
% 
%      stats=pyGAM.m3{i,1};
%     stats=struct(stats.statistics_);
%     stats.pseudo_r2=struct(stats.pseudo_r2);
%     stats.cov=double(stats.cov);
%     stats.p_values=struct(stats.p_values);
%     stats.edof_per_coef=struct(stats.edof_per_coef);
%     stats.m_features=int64(stats.m_features);
%     stats.n_samples=int64(stats.n_samples);
%     stats.se=double(stats.se);
%     structGAM.m3{i}=stats;
%     end
% end
% 
% cd '\\martinezsrv.robarts.ca\martinez_data$\Diego Buitrago-Piza\Data\Results\pyGAM'
% 
% save('onlyspeed_v2_struct_gamC.mat','structGAM')
% 
% save('gz_pl_hd_struct_gamPB.mat','structGAM')

%save('temp2_struct_gamC.mat','structGAM')