%% Read in the Ilastik prediction maps for membrane, convert h5 file to binnary tif file

close all;clear;clc;
for th=[0.2:0.1:0.9]; %probability threshold, you can specify one threshold that gives best segmentation
Amin = 100; Amax = 5000; %minimum and maximum area
ROI_name = '2-1_ROI02'; %change
dir_path = ['/Users/weichen/Desktop/Jan Lab/Kai_lyso_calcium/images/2-1 crop/',ROI_name,'/'];

im_name_base = [ROI_name,'_lys'];

h5name = [im_name_base,'_Probabilities.h5'];

ilastik_filename = [dir_path,h5name]; % ilastik prediction map here

im_lys = [dir_path,im_name_base,'.tif'];

info = imfinfo(im_lys);

Tmax=size(info,1);

stack_cal={};
stack_lys={};
for i=1:Tmax
imstack_lys{i}=imread(im_lys,i);
end

dir_write = [dir_path,im_name_base,'_segmentation/'];
delete(dir_write);
mkdir(dir_write);

ilastik_file = h5read(ilastik_filename,'/exported_data/'); % exported data is the "folder" in Ilastik project where the data is stored. Don't change this.
% Originally 4D file, we want 3D
label1 = squeeze(ilastik_file(1,:,:,:));
label2 = squeeze(ilastik_file(2,:,:,:)); 

label1 = permute(label1,[2,1,3]);                              % prevent the Ilastik axis swap by fixing them here
label1_db = double (label1);

label2 = permute(label2,[2,1,3]);                              % prevent the Ilastik axis swap by fixing them here
label2_db = double (label2);

seg = label1_db>th;
%% Apply area filters to bw image
    seg_new = seg;
for m = 1:Tmax
    seg_tmp = seg(:,:,m);
    seg_tmp = imfill(seg_tmp,'holes'); %fill holes inside objects
    tmp_lys = imstack_lys{m};
    stats= regionprops(seg_tmp,'all');
    obj_num = size(stats,1);

    for n = 1:obj_num
        if stats(n).Area<Amin || stats(n).Area > Amax
            Pxl_Id = stats(n).PixelIdxList;
            seg_tmp(Pxl_Id) = 0;
        else 
            continue;
        end
    end
    seg_new(:,:,m) = seg_tmp;
end

  %% save filtered tiff stack
im_write_name = [dir_write,'bw_th',num2str(th),'Amin',num2str(Amin),'Amax',num2str(Amax),'_',im_name_base,'.tif'];
    if exist(im_write_name)
        delete(im_write_name);
    end  
for ii = 1:Tmax;
    im0 = double(seg_new(:,:,ii));
    imwrite(im0,im_write_name,'tif','WriteMode','append');
end

end