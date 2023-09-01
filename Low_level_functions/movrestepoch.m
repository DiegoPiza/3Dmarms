function [mov rest start_end_mov  peaks_mov start_end_rest  peaks_rest locs amp_mov ]=movrestepoch(TrackingData)
%%
% This function processes tracking data to extract movement and resting epochs.
% The function returns various statistics related to the movement and resting data.

%movement epochs

b_vel=TracVelocity(TrackingData,1/60,0);

epoch_b=find_epochs(b_vel,16,true,12); %joins epochs whose start and end are within 200ms


amplitudeb=nan([length(epoch_b) 1]);
pksb=nan([length(epoch_b) 1]);
onsetb=nan([length(epoch_b) 1]);
termb=nan([length(epoch_b) 1]);
locsb=nan([length(epoch_b) 1]);
lengthb=nan([length(epoch_b) 1]);



for i=1:length(epoch_b)

% a=min_speed(find(min_speed<epoch_ang(i).start_time));
 onsetb(i)=epoch_b(i).start_time;
% b=min_speed(find(min_speed>epoch_ang(i).end_time));
if ~isempty(epoch_b(i).end_time)
    termb(i)=epoch_b(i).end_time;
else
   termb(i)=length(b_vel);
end
    [pksb(i) locsb(i)]=max(b_vel(onsetb(i): termb(i)));
    locsb(i)=onsetb(i)+locsb(i);
    amplitudeb(i)=sum((b_vel(onsetb(i):termb(i)).*0.0167));
if ~isempty(epoch_b(i).end_time)
    lengthb(i)=epoch_b(i).epoch_length;
else
    lengthb(i)=termb(i)-onsetb(i);
end
end


% 
% [h_vel v_vel] = HeadVelocity(TrackingData,sr,0) ;
% % view=Head_Orientation_No_floors(TrackingData);
% %south=0 north=1 , left=0, right=1
% 
% 
% Threshdiff = [2];
% PeakThresh = [3000];
%     for var = 1
%         
%         while Threshdiff(var)>1
%             MeanVel(var) = mean(b_vel(b_vel(:,var)<PeakThresh(var),var));
%             STD(var) = std(b_vel(b_vel(:,var)<PeakThresh(var),var));
%             NewThreshb = MeanVel(var)+(sigmab*STD(var));
%             Threshdiff(var) = PeakThresh(var)-NewThreshb;
%             PeakThresh(var) = NewThreshb;
%             i=i+1;
%             
%         end
%         
%     end
%     
% 
% min_speedb = islocalmin(b_vel);
% min_speedb=find(min_speedb);
% 
% NewThreshb=10;
% 
% [pksb,locsb] =findpeaks(b_vel,'MinPeakHeight',NewThreshb);
% 
% last=min_speedb(find(min_speedb>locsb(end)));
% if isempty(last)
%     pksb=pksb(1:end-1);
%     locsb=locsb(1:end-1);
% end
% 
% first= min_speedb(find(min_speedb<locsb(1)));
% if isempty(first)
%     pksb=pksb(2:end);
%     locsb=locsb(2:end);
% end
% %
% onsetb=nan([length(locsb) 1]);
% termb=nan([length(locsb) 1]);
% amplitudeb=nan([length(locsb) 1]);
% for i=1:length(locsb)
% a=min_speedb(find(min_speedb<locsb(i)&b_vel(min_speedb)<(pksb(i)/2)));
% if isempty(a)
%     a=min_speedb(find(min_speedb<locsb(i))); 
% end
% onsetb(i)=a(end);
% b=min_speedb(find(min_speedb>locsb(i)&b_vel(min_speedb)<(pksb(i)/2)));
% termb(i)=b(1);
% amplitudeb(i)=sum((b_vel(onsetb(i):termb(i)).*0.0167));
% % if pksb(i)>600
% % x=TrackingData.kZQuat(onsetb(i)-10:term(i)+10)  ;
% % if x>100
% %     continue
% % end
% % % x=rescale(x,0,1);
% % % plot(x,'b');
% % % hold on
% % end
% end
% 
onsetb(amplitudeb<15)=[];
termb(amplitudeb<15)=  [];                %movements less than 30cm in length
pksb(amplitudeb<15)=[];
locsb(amplitudeb<15)=[];
lengthb(amplitudeb<15)=[];
amplitudeb(amplitudeb<15)=[];


onsetb(lengthb<30)=[];
termb(lengthb<30)=  [];                %movements less than 30cm in length
pksb(lengthb<30)=[];
locsb(lengthb<30)=[];
amplitudeb(lengthb<30)=[];
lengthb(lengthb<30)=[];




% todelete=[];
% index=[];
% for i=1:length(onset)
% if sum(isnan(ang_vel(onset(i):term(i))))>0
%     todelete=[todelete;i];
%     continue
% else
% index=[index onset(i):term(i)];
% end
% end
% locs(todelete,:)=[];
% onset(todelete,:)=[];
% term(todelete,:)=[];
%
todelete=[];
indexb=[];
df=nan(size(onsetb));
df(2:end)=diff(onsetb);
for i=1:length(onsetb)
if df(i)>0   
if sum(isnan(b_vel(onsetb(i):termb(i))))>0
    todelete=[todelete;i];
    continue
else
indexb=[indexb onsetb(i):termb(i)];
pksb(i)=max(b_vel(onsetb(i):termb(i)));
end
else
    todelete=[todelete;i];
end
end
%
% index(todelete,:)=[];
locsb(todelete,:)=[];
onsetb(todelete,:)=[];
termb(todelete,:)=[];
pksb(todelete,:)=[];
amplitudeb(todelete,:)=[];

locs=locsb;

start_end_mov(1,:)=onsetb;
start_end_mov(2,:)=termb;
peaks_mov=pksb;
mov=indexb;
amp_mov=amplitudeb;




% rest epochs

epoch_b=find_epochs(b_vel,16,false,12); %joins epochs whose start and end are within 200ms


amplitudeb=nan([length(epoch_b) 1]);
pksb=nan([length(epoch_b) 1]);
onsetb=nan([length(epoch_b) 1]);
termb=nan([length(epoch_b) 1]);
locsb=nan([length(epoch_b) 1]);
lengthb=nan([length(epoch_b) 1]);



for i=1:length(epoch_b)

% a=min_speed(find(min_speed<epoch_ang(i).start_time));
 onsetb(i)=epoch_b(i).start_time;
% b=min_speed(find(min_speed>epoch_ang(i).end_time));
if ~isempty(epoch_b(i).end_time)
    termb(i)=epoch_b(i).end_time;
else
   termb(i)=length(b_vel);
end
    [pksb(i) locsb(i)]=max(b_vel(onsetb(i): termb(i)));
    amplitudeb(i)=sum((b_vel(onsetb(i):termb(i)).*0.0167));
if ~isempty(epoch_b(i).end_time)
    lengthb(i)=epoch_b(i).epoch_length;
else
    lengthb(i)=termb(i)-onsetb(i);
end
end


% 
% [h_vel v_vel] = HeadVelocity(TrackingData,sr,0) ;
% % view=Head_Orientation_No_floors(TrackingData);
% %south=0 north=1 , left=0, right=1
% 
% 
% Threshdiff = [2];
% PeakThresh = [3000];
%     for var = 1
%         
%         while Threshdiff(var)>1
%             MeanVel(var) = mean(b_vel(b_vel(:,var)<PeakThresh(var),var));
%             STD(var) = std(b_vel(b_vel(:,var)<PeakThresh(var),var));
%             NewThreshb = MeanVel(var)+(sigmab*STD(var));
%             Threshdiff(var) = PeakThresh(var)-NewThreshb;
%             PeakThresh(var) = NewThreshb;
%             i=i+1;
%             
%         end
%         
%     end
%     
% 
% min_speedb = islocalmin(b_vel);
% min_speedb=find(min_speedb);
% 
% NewThreshb=10;
% 
% [pksb,locsb] =findpeaks(b_vel,'MinPeakHeight',NewThreshb);
% 
% last=min_speedb(find(min_speedb>locsb(end)));
% if isempty(last)
%     pksb=pksb(1:end-1);
%     locsb=locsb(1:end-1);
% end
% 
% first= min_speedb(find(min_speedb<locsb(1)));
% if isempty(first)
%     pksb=pksb(2:end);
%     locsb=locsb(2:end);
% end
% %
% onsetb=nan([length(locsb) 1]);
% termb=nan([length(locsb) 1]);
% amplitudeb=nan([length(locsb) 1]);
% for i=1:length(locsb)
% a=min_speedb(find(min_speedb<locsb(i)&b_vel(min_speedb)<(pksb(i)/2)));
% if isempty(a)
%     a=min_speedb(find(min_speedb<locsb(i))); 
% end
% onsetb(i)=a(end);
% b=min_speedb(find(min_speedb>locsb(i)&b_vel(min_speedb)<(pksb(i)/2)));
% termb(i)=b(1);
% amplitudeb(i)=sum((b_vel(onsetb(i):termb(i)).*0.0167));
% % if pksb(i)>600
% % x=TrackingData.kZQuat(onsetb(i)-10:term(i)+10)  ;
% % if x>100
% %     continue
% % end
% % % x=rescale(x,0,1);
% % % plot(x,'b');
% % % hold on
% % end
% end
% 
% onsetb(amplitudeb<15)=[];
% termb(amplitudeb<15)=  [];                %movements less than 30cm in length
% pksb(amplitudeb<15)=[];
% locsb(amplitudeb<15)=[];
% lengthb(amplitudeb<15)=[];
% amplitudeb(amplitudeb<15)=[];


onsetb(lengthb<30)=[];
termb(lengthb<30)=  [];                %movements less than 30 samples
pksb(lengthb<30)=[];
locsb(lengthb<30)=[];
amplitudeb(lengthb<30)=[];
lengthb(lengthb<30)=[];


% todelete=[];
% index=[];
% for i=1:length(onset)
% if sum(isnan(ang_vel(onset(i):term(i))))>0
%     todelete=[todelete;i];
%     continue
% else
% index=[index onset(i):term(i)];
% end
% end
% locs(todelete,:)=[];
% onset(todelete,:)=[];
% term(todelete,:)=[];
%
todelete=[];
indexb=[];
df=nan(size(onsetb));
df(2:end)=diff(onsetb);
for i=1:length(onsetb)
if df(i)>0   
if sum(isnan(b_vel(onsetb(i):termb(i))))>0
    todelete=[todelete;i];
    continue
else
indexb=[indexb onsetb(i):termb(i)];
pksb(i)=max(b_vel(onsetb(i):termb(i)));
end
else
    todelete=[todelete;i];
end
end
%
% index(todelete,:)=[];
locsb(todelete,:)=[];
onsetb(todelete,:)=[];
termb(todelete,:)=[];
pksb(todelete,:)=[];
amplitudeb(todelete,:)=[];

rest=indexb;


% deleting repeat values

a=ismember(rest,mov);
rest(a)=[];


start_end_rest(1,:)=onsetb;
start_end_rest(2,:)=termb;
peaks_rest=pksb;


end
