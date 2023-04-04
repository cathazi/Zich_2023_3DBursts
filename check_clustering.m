function check_clustering(face_center,SS_ROI_ind,cl_ft,cl_surf,cl_surf_org,t_ind,f_ind,e,c)

% delete png images of previous clusters
PATH_CL = '/Users/Catha/Desktop/UCL_P2_45/vis_cl/';
cd(PATH_CL);delete *.png
%cols = [31,120,180;51,160,44;227,26,28;255,127,0;202,178,214;106,61,154;255,255,153;177,89,40;166,206,227;178,223,138;251,154,153;253,191,111]./255;

after_tf_cl = unique(cl_surf);after_tf_cl(after_tf_cl==0)=[];
after_cols = hsv(length(after_tf_cl));
for z = 1:length(t_ind)
    f2 = figure('units','normalized','outerposition',[0 0 0.65 0.65],'DefaultAxesFontSize',18);hold on;
    
    % before tf clustering
    subplot(1,5,[1 2]);hold on;
    data = cl_surf_org(:,f_ind(z),t_ind(z));
    bef_tf_cl = unique(data);bef_tf_cl(bef_tf_cl==0)=[];
    bef_cols = flipud(hsv(length(bef_tf_cl)));
    for u = 1:length(bef_tf_cl)
        mask = data == bef_tf_cl(u);
        scatter3(face_center(SS_ROI_ind(mask),1),face_center(SS_ROI_ind(mask),2),face_center(SS_ROI_ind(mask),3),36,bef_cols(u,:),'filled');
    end
    grid on;axis square;view(-39,27);
    xlim([-80 80]);ylim([-100 100]);zlim([-60 100]);legend;
    title('before tf clustering');
    
    % after tf clustering
    subplot(1,5,[3 4]);hold on;
    col_counter = [];
    for u = 1:length(after_tf_cl)
        data = cl_surf(:,f_ind(z),t_ind(z));
        mask = data == after_tf_cl(u);
        if sum(mask)>0
            scatter3(face_center(SS_ROI_ind(mask),1),face_center(SS_ROI_ind(mask),2),face_center(SS_ROI_ind(mask),3),36,after_cols(u,:),'filled');
            col_counter{end+1} = num2str(u);
        end
    end
    grid on;axis square;view(-39,27);
    xlim([-80 80]);ylim([-100 100]);zlim([-60 100]);legend;legend(col_counter);
    title('after tf clustering');
    
    % tf matrix
    h=subplot(1,5,5);hold on;
    imagesc(cl_ft);colormap(h,flip(hot));caxis([0 1]);
    plot([t_ind(z) t_ind(z)],[1 18],'m','Linewidth',2);
    plot([1 72],[f_ind(z) f_ind(z)],'c','Linewidth',2);
    axis square;xlim([1 72]);ylim([1 18]);
    xticks([1:10:72]);yticks([1:2:18]);yticklabels([13:2:30]);xlabel('Time (samples)');ylabel('Frequency (Hz)');
    title(['Freq: ',num2str(12+f_ind(z)),'/ Time: ',num2str(t_ind(z))]);
    
    print('-dpng',[PATH_CL 'Epoch_' num2str(e,'%03.f') '_cl_' num2str(c,'%02.f') '_t_' num2str(t_ind(z),'%02.f') '_f_' num2str(f_ind(z),'%02.f') '.png']);
    close(f2);
end


get_surf_cluster_check_vid(PATH_CL,e,c)