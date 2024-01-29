function [fix, start_end_sac, peaks_sac, start_end_fix, peaks_fix, locss, amplitudesac] = fixsacepoch(TrackingData)
%FIXSACEPOCH Process tracking data to determine rapid head movements and fixation epochs
%  fix= fixation timestamps
% start_end_sac= onset and offset timestamps of head movements
% peaks_sac= peak head movement velocity 
% start_end_fix= onset and offset timestamps of head fixations
% locss= timestamps for the peak velocity of head movements
% amplitudesac= amplitude in degrees of head movement

    % Calculate angular velocity
ang_vel=HeadSpeed(TrackingData,1/60,0);

[h_vel v_vel r_vel] = HeadVelocity(TrackingData,1/60,0);


sr=60; %sampling rate

    % Find saccade epochs

epoch_ang=find_epochs(ang_vel,200,true,4);
    % Initialize variables

amplitude=nan([length(epoch_ang) 1]);
pks=nan([length(epoch_ang) 1]);
onset=nan([length(epoch_ang) 1]);
term=nan([length(epoch_ang) 1]);
locs=nan([length(epoch_ang) 1]);
lengtha=nan([length(epoch_ang) 1]);

    % Extract parameters for saccade epochs

for i=1:length(epoch_ang)

    onset(i)=epoch_ang(i).start_time;
    if ~isempty(epoch_ang(i).end_time)
        term(i)=epoch_ang(i).end_time;
    else
        term(i)=length(ang_vel);
    end
    [pks(i) locs(i)]=max(ang_vel(onset(i): term(i)));
    locs(i)=onset(i)+locs(i);
    amplitude(i)=sum((ang_vel(onset(i):term(i)).*0.0167));
    if ~isempty(epoch_ang(i).end_time)
        lengtha(i)=epoch_ang(i).epoch_length;
    else
        lengtha(i)=term(i)-onset(i);
    end
end

    % Filter out movements with less than 10 deg amplitude and invalida frames

onset(amplitude<10)=[];
term(amplitude<10)=  [];                %movements less than 10 deg amplitude
pks(amplitude<10)=[];
locs(amplitude<10)=[];
lengtha(amplitude<10)=[];
amplitude(amplitude<10)=[];



todelete=[];
index=[];
for i=1:length(onset)
    if sum(isnan(ang_vel(onset(i):term(i))))>0
        todelete=[todelete;i];
        continue
    else
        index=[index onset(i):term(i)];
    end
end
locs(todelete,:)=[];
onset(todelete,:)=[];
term(todelete,:)=[];
pks(todelete,:)=[];
amplitude(todelete,:)=[];
todelete=[];
index=[];
for i=1:length(onset)
    if sum(isnan(ang_vel(onset(i):term(i))))>0
        todelete=[todelete;i];
        continue
    else
        index=[index onset(i):term(i)];
    end
end
locs(todelete,:)=[];
onset(todelete,:)=[];
term(todelete,:)=[];
pks(todelete,:)=[];
amplitude(todelete,:)=[];
amplitudesac=amplitude;

  

sac=index;


start_end_sac(1,:)=onset;
start_end_sac(2,:)=term;
peaks_sac=pks;
locss=locs;

directionv=nan(size(onset));
directionh=nan(size(onset));

directionv(v_vel(locs)>0)=1; %Up
directionv(v_vel(locs)<0)=-1; % Down
directionh(h_vel(locs)>0)=1; %left
directionh(h_vel(locs)<0)=-1; % right

%fix

  % Find fixation epochs
epoch_ang=find_epochs(ang_vel,200,false,1);
amplitude=nan([length(epoch_ang) 1]);
pks=nan([length(epoch_ang) 1]);
onset=nan([length(epoch_ang) 1]);
term=nan([length(epoch_ang) 1]);
locs=nan([length(epoch_ang) 1]);
lengtha=nan([length(epoch_ang) 1]);
    % Extract parameters for fixation epochs

for i=1:length(epoch_ang)

    onset(i)=epoch_ang(i).start_time;
    if ~isempty(epoch_ang(i).end_time)
        term(i)=epoch_ang(i).end_time;
    else
        term(i)=length(ang_vel);
    end
    [pks(i) locs(i)]=max(ang_vel(onset(i): term(i)));
    amplitude(i)=sum((ang_vel(onset(i):term(i)).*0.0167));
    if ~isempty(epoch_ang(i).end_time)
        lengtha(i)=epoch_ang(i).epoch_length;
    else
        lengtha(i)=term(i)-onset(i);
    end
end

    % Filter out fixations with less than 200ms in duration

time_cutoff=floor(sr*0.2);
onset(lengtha<time_cutoff)=[];
term(lengtha<time_cutoff)=  [];                %fix less than 200ms in duration
pks(lengtha<time_cutoff)=[];
locs(lengtha<time_cutoff)=[];
amplitude(lengtha<time_cutoff)=[];
lengtha(lengtha<time_cutoff)=[];


todelete=[];
index=[];
for i=1:length(onset)
    if sum(isnan(ang_vel(onset(i):term(i))))>0
        todelete=[todelete;i];
        continue
    else
        if   length(onset(i)+1:term(i)-1)<11
            error('not correct length')
        end
        index=[index onset(i)+1:term(i)-1];
    end
end
locs(todelete,:)=[];
onset(todelete,:)=[];
term(todelete,:)=[];
pks(todelete,:)=[];



fix=index;

a=ismember(fix,sac);
fix(a)=[];


start_end_fix(1,:)=onset;
start_end_fix(2,:)=term;
peaks_fix=pks;
end
