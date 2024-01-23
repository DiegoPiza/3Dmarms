%% Speed Time series analysis  https://www.sciencedirect.com/science/article/pii/S2211124718316437#fig3
%Zé Henrique T.D. Góis - Characterizing Speed Cells in the Rat Hippocampus
rng shuffle
clear all
close all
%%
load sample_data.mat
%%
speed_score=table;

numofshuffle=1000; %number of shuffle
sr=60; %sampling rate of  time series

%gaussian kernel
sigma=floor(0.25*sr);%250 ms window
L=sigma*12; %length of 6 sigma each side
alpha=((L-1)/sigma)/2;
h = 1/(sqrt(2*pi)*sigma);   % height of Gaussian
w=gausswin(L,alpha)*h; % The exact correspondence with the standard deviation of a Gaussian probability density function is σ = (L – 1)/(2α).
w(w < h*1.e-3)=[]; %deleting gausian edges
plot_figures=false;

fields=fieldnames(units.singleunits);
score=table;      
for j=1:length(fields)
unitID=string(fields(j));
neuron=char(fields(j))

unit=units.singleunits.(neuron);
if length(unit)<250
    continue
end
score.unitID{j}=unitID;

%%
puls=TrackingData.neural_time;
close all

gapt=find(gap)./30000; %in seconds

inv_frm=ismember(puls,gapt);

spikes=histc(unit,puls);
heatm=HeadSpeed(TrackingData,1/60,1);
heatm(TrackingData.XPosition>0.6,:)=nan;
heatm(TrackingData.XPosition<0,:)=nan;
heatm(TrackingData.YPosition<0,:)=nan; 
heatm(TrackingData.YPosition>1.14,:)=nan;
heatm(TrackingData.ZPosition<0,:)=nan;

heatm(heatm>2000,:)=nan;

if length(spikes)> length(heatm) %Cheks if Tracking data is longer than the neural data
spikes=spikes(1:length(heatm));
inv_frm=inv_frm(1:length(heatm));
else 
heatm=heatm(1:length(spikes),:); 
inv_frm=inv_frm(1:length(spikes),:);
end

heatm(inv_frm,:)=nan;

%gausian kernel
spikes=conv(spikes,w,'same'); 
%deleting nans
spikes(isnan(heatm(:,1)))=[];
heatm(isnan(heatm(:,1)),:)=[];


R = corr(spikes,heatm)
score.Head_real{j}=R;

Rsh=zeros([numofshuffle 1]);
shift=randi([-90*sr 90*sr],[1000,1]); %-+ 90 sec shift

parfor i=1:numofshuffle
ash=circshift(spikes,shift(i)); 
Rsh(i)=corr(ash,heatm);
end
score.Head_sh{j}=Rsh;
%statistics
    mu = mean(Rsh);
    sigma_stats = std(Rsh);
    p_val = 1 - normcdf(R,mu,sigma_stats);
    reject_h0 = p_val < 0.05/length(spikes);
score.head_p_val{j}=p_val;
score.head_reject_h0{j}=reject_h0;       
%% plotting

if plot_figures==true
            h_fig2 = figure('units','normalized','position',[0 0 1 1]);
            subplot(2,1,1)
            histogram(Rsh,'normalization','pdf');
             hold all
             shuff_xs = linspace(min(Rsh),max(Rsh),100);
             shuff_curve = normpdf(shuff_xs,mu,sigma_stats);
             h_line1 = plot(shuff_xs,shuff_curve,'--b','linewidth',2);
             h_line2 = line([R R],get(gca,'ylim'),'color','r','linewidth',3);
             legend([h_line1,h_line2],{'Shuffled'; 'Real'},'location','northeastoutside')
             title('Head Speed Distribution of shufled total Speed Score')
             xlabel({'Speed Score (Corr)'})
             ylabel('Probability density')
             text(0.2,10,char(['Head Speed Score: ' num2str(R)]))
            text(0.2,15,char(['Permutation test p-value: ' num2str(p_val)]))
             
end



spikes=histc(unit,puls);

heatm=TracVelocity(TrackingData,1/60,1);
heatm(TrackingData.XPosition>0.6,:)=nan;
heatm(TrackingData.XPosition<0,:)=nan;
heatm(TrackingData.YPosition<0,:)=nan; 
heatm(TrackingData.YPosition>1.14,:)=nan;
heatm(TrackingData.ZPosition<0,:)=nan;
heatm(heatm>300,:)=nan;


if length(spikes)> length(heatm) %Cheks if Tracking data is longer than the neural data
spikes=spikes(1:length(heatm));
inv_frm=inv_frm(1:length(heatm));
else 
heatm=heatm(1:length(spikes),:); 
inv_frm=inv_frm(1:length(spikes),:);
end

heatm(inv_frm,:)=nan;



%gausian kernel
spikes=conv(spikes,w,'same'); 
%deleting nans
spikes(isnan(heatm(:,1)))=[];
heatm(isnan(heatm(:,1)),:)=[];

R = corr(spikes,heatm);
score.Body_real{j}=R;

Rsh=zeros([numofshuffle 1]);
shift=randi([-90*sr 90*sr],[1000,1]); %-+ 90 sec shift

parfor i=1:numofshuffle
ash=circshift(spikes,shift(i)); 
Rsh(i)=corr(ash,heatm);
end


score.Body_sh{j}=Rsh;
%statistics
    mu = mean(Rsh);
    sigma_stats = std(Rsh);
    p_val = 1 - normcdf(R,mu,sigma_stats);
    reject_h0 = p_val < 0.05/length(spikes);
score.body_p_val{j}=p_val;
score.body_reject_h0{j}=reject_h0;
         


%% plotting

if plot_figures==true
            subplot(2,1,2)
            
            histogram(Rsh,'normalization','pdf');
             hold all
             shuff_xs = linspace(min(Rsh),max(Rsh),100);
             shuff_curve = normpdf(shuff_xs,mu,sigma_stats);
             h_line1 = plot(shuff_xs,shuff_curve,'--b','linewidth',2);
             h_line2 = line([R R],get(gca,'ylim'),'color','r','linewidth',3);
             legend([h_line1,h_line2],{'Shuffled'; 'Real'},'location','northeastoutside')
             title('Body Speed Distribution of shufled total Speed Score')
             xlabel({'Speed Score (Corr)'})
             ylabel('Probability density')
             text(0.2,10,char(['Body Vel Speed Score: ' num2str(R)]))
            text(0.2,15,char(['Permutation test p-value: ' num2str(p_val)]))           
                       
end

end

