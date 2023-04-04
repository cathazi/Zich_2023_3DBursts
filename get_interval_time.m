% get interval time

function inti = get_interval_time(cl_info)

% get index of Epochs containing clusters
Epo = unique(cl_info.epoch_index);

counter = 0;
for e = 1:length(Epo)
    
    % get on- offset of clusters
    cl_inEpo = find(cl_info.epoch_index==Epo(e));
    clear onoff
    onoff(:,1) =  cl_info.t_onset(cl_inEpo);
    onoff(:,2) =  cl_info.t_offset(cl_inEpo);
    
    onoff(onoff(:,1)==0,:)=[];
    
    % sort clusters based on onset time
    onoff = sortrows(onoff,1);
    for i = 1:size(onoff,1)-1
        counter = counter+1;
        inti(counter) = onoff(i+1,1) - onoff(i,2);
    end
end % Epochs
