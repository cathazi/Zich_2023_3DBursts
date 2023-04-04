% creates video out of images

function get_surf_cluster_check_vid(PATH_CL,e,c)

cd(PATH_CL);
imgs = dir('*.png');
writerObj = VideoWriter(['Epoch_' num2str(e,'%03.f') '_cl_' num2str(c,'%02.f') '.avi']);
writerObj.FrameRate = 2;
open(writerObj);
for i = 1 : length(imgs)
    thisimage = imread(imgs(i).name);
    writeVideo(writerObj, thisimage);
end
close(writerObj);