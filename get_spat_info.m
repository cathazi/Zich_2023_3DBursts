function [sl,sl_vi,SS_ROI_ind,ROI_ind] = get_spat_info(PATH_SUBJ,PATH_SURF,PATH_MRI,SUBJ)

% load surface
cd(PATH_SURF);surface = dir('S*mesh.surf.gii');temp = gifti([surface.folder '/' surface.name]);
% from gifti to struct
sl.faces = temp.faces;sl.vertices = temp.vertices;sl.mat = temp.mat;

surface = dir('S*mesh_vi.surf.gii');temp = gifti([surface.folder '/' surface.name]);
sl_vi.faces = temp.faces;sl_vi.vertices = temp.vertices;sl_vi.mat = temp.mat;

% get center of faces
for f = 1:size(sl.faces,1)
    sl.face_center(f,:) = double(mean(sl.vertices(sl.faces(f,:),:),1));
    sl_vi.face_center(f,:) = double(mean(sl_vi.vertices(sl_vi.faces(f,:),:),1));
end

% % check center of faces
% temp = ones(size(sl.faces,1),1);
% h=figure;hold on;
% p = patch('Faces',sl.faces,'Vertices',sl.vertices,'FaceVertexCData',temp,'EdgeColor',[1 1 1],'FaceColor','flat');
% scatter3(face_center(:,1),face_center(:,2),face_center(:,3),'r*');

% % check distance between center of faces
% [idx,dist] = knnsearch(sl.vertices,sl.vertices);
% figure;
% histogram(dist);xlabel('Distance (mm)');ylabel('N of occurances');title('distance between face centeres');

% load search space ROI
% Search space = frontal and parietal lobe (whole left hemisphere is too big to save after beamform, if needed, beaform each trial individually)
cd(PATH_MRI);MRI = dir('p*HeadCast*');PATH_MRI_SUBJ = [PATH_MRI MRI.name '/'];cd(PATH_MRI_SUBJ);SS_roi = dir('*SS*coords.mat');
load([PATH_MRI_SUBJ SS_roi.name]);
in = inhull(sl.face_center,roi_coords.mni',[],3);


% to check tolerance parameter in inhull
% figure;
% subplot(1,3,1);h=scatter3(roi_coords.mni(1,:),roi_coords.mni(2,:),roi_coords.mni(3,:),'k*');title('mask');xlim([-100 100]);ylim([-100 100]);zlim([-100 100]);axis square;
% subplot(1,3,2);h=scatter3(sl.vertices(in,1),sl.vertices(in,2),sl.vertices(in,3),'k*');title('captured');xlim([-100 100]);ylim([-100 100]);zlim([-100 100]);axis square;
% subplot(1,3,3);h=scatter3(sl.vertices(not(in),1),sl.vertices(not(in),2),sl.vertices(not(in),3),'k*');title('not captured');xlim([-100 100]);ylim([-100 100]);zlim([-100 100]);axis square;

SS_ROI_ind = in.*[1:length(in)]';SS_ROI_ind(SS_ROI_ind==0)=[];
% make sure they are not too close - otherwise beamformer problem
% currently remove too close once (results in holes in brain)
if strcmp(SUBJ,'rh') || strcmp(SUBJ,'rk') || strcmp(SUBJ,'ad_no_helmet')
    [idx,d] = knnsearch(sl.face_center(SS_ROI_ind,:),sl.face_center(SS_ROI_ind,:),'K',2);
    SS_ROI_ind(idx(d(:,2)<0.1)) = [];
end

% get indecies of ROI relative to SS_ROI (quality check of SS_ROI and ROI)
ROI_ind = get_ROI_ind(PATH_MRI_SUBJ,SUBJ,SS_ROI_ind,sl,sl_vi,0);

% save spatial params per subj
save([PATH_SUBJ 'spat_info.mat'],'sl','sl_vi','SS_ROI_ind','ROI_ind','-v7.3');
