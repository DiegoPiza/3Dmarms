function [mov rest start_end_mov  peaks_mov start_end_rest  peaks_rest locs amp_mov ]=movrestepoch(TrackingData)
%%
% This function processes tracking data (assumes sampling rate 60 Hz) to extract movement and resting epochs.
% The function returns various statistics related to the movement and resting data.
%  mov= movement timestamps
% rest= no movement timestamps
% start_end_mov= onset and offset timestamps of movements
% peaks_mov= peak movement velocity 
% start_end_rest= onset and offset timestamps of rest epochs
% peaks_rest= peak speed of  rest epochs
% locs= timestamps for the peak velocity of head movements
% amp_mov= amplitude in cm movements

%movement epochs

b_vel=TracVelocity(TrackingData,1/60,0);

epoch_b=find_epochs(b_vel,15,true,12); %joins epochs whose start and end are within 200ms


amplitudeb=nan([length(epoch_b) 1]);
pksb=nan([length(epoch_b) 1]);
onsetb=nan([length(epoch_b) 1]);
termb=nan([length(epoch_b) 1]);
locsb=nan([length(epoch_b) 1]);
lengthb=nan([length(epoch_b) 1]);



for i=1:length(epoch_b)

 onsetb(i)=epoch_b(i).start_time;
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


onsetb(amplitudeb<15)=[];
termb(amplitudeb<15)=  [];                %movements less than 15cm in length
pksb(amplitudeb<15)=[];
locsb(amplitudeb<15)=[];
lengthb(amplitudeb<15)=[];
amplitudeb(amplitudeb<15)=[];


onsetb(lengthb<30)=[];
termb(lengthb<30)=  [];                %movements less than 0.5s duration
pksb(lengthb<30)=[];
locsb(lengthb<30)=[];
amplitudeb(lengthb<30)=[];
lengthb(lengthb<30)=[];

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

epoch_b=find_epochs(b_vel,15,false,12); %joins epochs whose start and end are within 200ms

amplitudeb=nan([length(epoch_b) 1]);
pksb=nan([length(epoch_b) 1]);
onsetb=nan([length(epoch_b) 1]);
termb=nan([length(epoch_b) 1]);
locsb=nan([length(epoch_b) 1]);
lengthb=nan([length(epoch_b) 1]);

for i=1:length(epoch_b)

 onsetb(i)=epoch_b(i).start_time;
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

onsetb(lengthb<30)=[];
termb(lengthb<30)=  [];                %movements less than 30 samples or 0.5s duration
pksb(lengthb<30)=[];
locsb(lengthb<30)=[];
amplitudeb(lengthb<30)=[];
lengthb(lengthb<30)=[];
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
