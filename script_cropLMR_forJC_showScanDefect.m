% chenzhe, 2017-01-18
% crop the 1st, 4th, 5th, and 8th stack
% regions on the left, middle, and right
% to illustrate the scan defect
% The cropped image will be given to JC-stinville

oldDir = pwd;

targetDir = 'D:\UMich Folder\11 ONR project\32 Microscopy\2017-11-07 External illustrate distortion\imgs_s4_n1_l1_snake\individual_stack_s4_n1_l1_snake';
cd(targetDir);

cropRegion = {[0, 1200, 48, 48];
    [940, 1200, 48, 48];
    [2000, 1200, 48, 48]};

for iStack=[1,4,5,8]
    imgName = strcat('frame_0',num2str(iStack),'.tif');
    img = imread(imgName);
    for iArea=1:3
        imgCropped = imcrop(img,cropRegion{iArea});
        % change the stack number according to the physical scan sequence
        if iStack==5
            iStackM = 8;
        elseif iStack==8
            iStackM = 5;
        else
            iStackM = iStack;
        end
        outputNname = strcat('area_',num2str(iArea),'_stack_',num2str(iStackM),'.tif');
        imwrite(imgCropped,outputNname);
    end
end
cd(oldDir);