function data = check_sign_ambiguity(data,c,vis)
% add: show data before and after sign flip?

% visualise time course of each voxel
if vis
    f1 = figure('DefaultAxesFontSize',16);
    subplot(2,3,[1 2]);hold on;
    plot(data');
    plot(mean(data,1)','k','Linewidth',3);
    plot([1 size(data,2)],[0 0],'r');
    axis tight;
    xlabel('time (samples)');ylabel('power (db)');
    title(['Cluster: ',num2str(c),', N of voxels: ',num2str(size(data,1))]);
end

% get correlation with avg. time series (avg. might not be ideal if 50:50 split)
for v = 1:size(data,1)
    r_temp = corrcoef(mean(data,1),data(v,:));
    r(v) = r_temp(1,2);
end

% visualise correlation coefficients
if vis
    subplot(2,3,3);hold on;
    plot(r,'k');
    xlabel('voxels');ylabel('correalation coef');
    
end

% flip if necessary
if sum(r<0)>0
    data((r<0),:) = data((r<0),:)*-1;
    
    if vis
        % visualise time course of each voxel (post flip)
        subplot(2,3,[4 5]);hold on;
        plot(data');
        plot(mean(data,1)','k','Linewidth',3);
        plot([1 size(data,2)],[0 0],'r');
        axis tight;
        xlabel('time (samples)');ylabel('power (db)');
        title('post sign flip');
    end
end

