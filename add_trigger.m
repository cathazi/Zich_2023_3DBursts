function [D_cont,eve_info,RT] = add_trigger(D_cont,PATH_SESS)


Events=D_cont.events;
% delete evebnts added in previous iterations of this script
to_be_deleted = [];
for e=1:length(Events)
    if Events(e).value==66
        to_be_deleted(end+1)=e;
    end
end
Events(to_be_deleted)=[];

% load response & stimuli information
parts = strsplit(D_cont.fname,{'_','.'});
load([PATH_SESS '/data_' parts{4} '_' num2str(str2num(parts{8}))]);
load([PATH_SESS '/stim_' parts{4} '_' num2str(str2num(parts{8}))]);

% adjust for the fact that 1/2 = L/R in data.responses to not match
% 1/2 = L/R in stim.trials(:,1) & stim.trials(:,4)
stim.trials(:,1)=1+2-stim.trials(:,1);
stim.trials(:,4)=1+2-stim.trials(:,4);

% get levels (different across individuals)
levels=unique(stim.trials(:,2));

% bring all relevant information in one variable
counter=0;
for e=1:size(data.responses,1)
    if data.responses(e,2)>0 % exclude trials with no response
        counter = counter+1;
        % RT (similar to the time between trigger 50 and 60 in MEG)
        eve_info(counter,1) = data.responses(e,2); 
        % cong/incong = 1/0
        eve_info(counter,2) = stim.trials(e,3); 
        % add level information (not real values, but 1/2/3 = low/medium/high)
        for l=1:length(levels)
            if stim.trials(e,2) == levels(l)
                eve_info(counter,3) = l;
            end
        end
        % add accuracy information (accurate/ not accurate = 1/0)
        if stim.trials(e,4)==data.responses(e,1)
            eve_info(counter,4) = 1;
        else
            eve_info(counter,4) = 0;
        end
    end
end

Events_new_start = struct([]);RT=[];
st50_ind=0;
for e=1:length(Events)
    if strcmp(Events(e).type,'UPPT001_up') && Events(e).value==50
        st50_ind(end+1) = st50_ind(end)+1;
        st50_t(st50_ind(end)) = Events(e).time;
    elseif strcmp(Events(e).type,'UPPT001_up') && Events(e).value==60
        Events_new_start(end+1).type   = 'RT';
        Events_new_start(end).value    = 66;
        Events_new_start(end).time     = Events(e).time;
        Events_new_start(end).duration = [];
        Events_new_start(end).offset   = 0;
        st60_t(st50_ind(end)) = Events(e).time;
    end
end
st50_ind(st50_ind==0)=[]; %remove first entry in vector (artifically introduced)

% remove trials without response 
no_RT2=find(st60_t==0);
st50_t(no_RT2)=[];st60_t(no_RT2)=[];
RT = st60_t-st50_t; % RT
RT(1)=[]; % remove first trial (as in Bonaiuto scripts, 181 entries, but only 180 trials)
Events_new_start(1)=[];

% check
if length(RT) ~= size(eve_info,1)
    disp('dimension mismatch between behavioural data aquired within the meg and in parallel to the !')
    f1=figure;plot(RT,'k');hold on;plot(eve_info(:,1),'r--');legend('RT','eve info');
    %decision=input('Enter 1 to remove additional trials in eve info strcuture: ');
    decision=1;
    close(f1);
    if decision == 1
        % if MEG terminated before experiemnt finished, it might be that
        % there are no 2 seconds follwoing the last recorded RT, thus
        % exclude the last recorded RT to ensure match between behav and
        % meg data
        RT = RT(1:length(RT)-1);
        eve_info=eve_info(1:length(RT),:);
        Events_new_start=Events_new_start(1:length(RT));
    end
end

% select events
for r = 1:size(eve_info,1) 
    if eve_info(r,2) == 1 && eve_info(r,3) == 3 && eve_info(r,4) == 1 
        keep(r) = r;
    end
end
keep(keep==0)=[];
Events_new_start = Events_new_start(keep);
RT = RT(keep);    

if ~isempty(Events_new_start)
    Events = [Events(:); Events_new_start(:)];
end
D_cont = D_cont.events(1,Events);
D_cont.save