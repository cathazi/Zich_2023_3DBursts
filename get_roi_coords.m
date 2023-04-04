% does: 
% 1) segments T1 to get deformation fiels
% 2) normalises AAL MNI ROIs into single subject space
% 3) extracts and saves ROI coords

% needs: 
% define search space ROI (SS ROI) and sensorimotor ROI (SM ROI) in WFU PickAtlas (save as nii file)
% spm toolbox marsbar
% funcitons extract_voxel_values & vox2mni 

function get_roi_coords(PATH_SCRIPT,PATH_MRI_SUBJ,PATH_ROI)

cd(PATH_MRI_SUBJ);def_field = dir('iy*.nii');
if length(def_field)==0
    % Segment T1 and save inverse deformation field
    load([PATH_SCRIPT 'spm_segment.mat']);
    cd(PATH_MRI_SUBJ);T1 = dir('s2*.nii');
    matlabbatch{1}.spm.spatial.preproc.channel.vols = {[PATH_MRI_SUBJ T1.name ',1']};
    spm_jobman('run',matlabbatch)
elseif length(def_field)>1
    disp('Error: more than one inverse deforamtion field exists.')
end

% Normalise AAL_MNI_ROIs to subject space using the inverse
% deformation field
% delete previous w* files in ROI folder 
cd(PATH_ROI);w_rois = dir('w*.nii');
if length(w_rois)>0
    for w = 1: length(w_rois)
        delete(strtrim([w_rois.folder '/' w_rois.name]));
    end
end

load([PATH_SCRIPT 'spm_normalise.mat']);
cd(PATH_MRI_SUBJ);def_field = dir('iy*.nii');
matlabbatch{1}.spm.spatial.normalise.write.subj.def = {[PATH_MRI_SUBJ def_field.name]};
cd(PATH_ROI);roi_org_nii = dir('*.nii');
for r = 1:length(roi_org_nii)
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample(r,:) = {[roi_org_nii(r).folder '/' roi_org_nii(r).name ',1']};
end
spm_jobman('run',matlabbatch)


cd(PATH_ROI);roi_w_nii = dir('w*.nii');
for r = 1:length(roi_w_nii)
    % move normalised ROIs to subj MRI folder
    movefile([roi_w_nii(r).folder '/' roi_w_nii(r).name],PATH_MRI_SUBJ);
    
    % transform ROIs form .nii to _roi.mat (marsbar objects)
    temp = strsplit(roi_w_nii(r).name,'.');
    mars_img2rois([PATH_MRI_SUBJ roi_w_nii(r).name], PATH_MRI_SUBJ, temp{1}, 'x');

    % get coords of roi
    [~,roi_coords,~] = extract_voxel_values([PATH_MRI_SUBJ temp{1} '_1_roi.mat'],[PATH_MRI_SUBJ temp{1} '.nii']);
    % #### sometimes _1_ sometimes only _ (figure out and fix)
    roi_coords = roi_coords.I;
    save([PATH_MRI_SUBJ temp{1} '_coords.mat'],'roi_coords');
end
