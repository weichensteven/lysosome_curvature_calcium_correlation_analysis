
clear; close all; clc;
Start = tic;
addpath('/Users/weichen/Desktop/Jan Lab/Kai_lyso_calcium/github_repo'); 
% file path for curvature measure and visualization function written by Preetham Manjunatha
% https://www.mathworks.com/matlabcentral/fileexchange/96982-curvature-measure-and-visualization
%% Inputs
im_name_base = '2-1_ROI02';
Amin = 100; Amax = 5000; %minimum and maximum area when converted to bw image
th=0.8; % Probability threshold
dir_root = ['/Users/weichen/Desktop/Jan Lab/Kai_lyso_calcium/images/2-1 crop/',im_name_base,'/']; %input folder path
dir_write1 =[dir_root,'Curv_Color_ Outline2/']; 
dir_write2 =[dir_root,'IntensityRatio_Color_Outline2/'];
%dir_write3 =[dir_root,'IntensityCal_Color_outline/']; %only calcium
%channel
dir_write4 =[dir_root,'IntensityROI_outline/'];
thickness = 2; % thickness of the ROI band for intensity measure
%
if isfolder(dir_write1)
rmdir(dir_write1,'s');
end
mkdir(dir_write1);

if isfolder(dir_write2)
rmdir(dir_write2,'s');
end
mkdir(dir_write2);
%}
%{
if isfolder(dir_write3)
rmdir(dir_write3,'s');
end
mkdir(dir_write3);
%}
%
if isfolder(dir_write4)
rmdir(dir_write4,'s');
end
mkdir(dir_write4);
%}

im_path_lys = [dir_root,im_name_base,'_lys.tif']; %file path for lysosome channel image
im_path_cal = [dir_root,im_name_base,'_cal.tif']; %file path for calcium channel image
im_path_bw = [dir_root,im_name_base,'_lys_segmentation/bw_th',num2str(th),'Amin',num2str(Amin),'Amax',num2str(Amax),'_',im_name_base,'_lys.tif']; 
%file path for binary image

export_resolution = 200;  % resolution of exported images
color_width = 3; % width of the color coded band
info = imfinfo(im_path_lys);
im_height = info(1).Height;
im_width = info(1).Width;
Tmax=size(info,1);
stack_cal={};
stack_lys={};
stack_bw={};  % Binary image stack of segmentation
stack_bwlabel={}; % object-labeled image stack
R_all = [];
for i=1:Tmax
stack_cal{i}=imread(im_path_cal,i);
stack_lys{i}=imread(im_path_lys,i);
stack_bw{i}=imread(im_path_bw,i);
stack_bwlabel{i}=bwlabel(stack_bw{i});
end

C_all = []; %curvature
C2_all = []; %intensity ratio
all_obj_area = [];
all_obj_lys_intensity = [];

for k = 1:Tmax
    im_tmp = stack_bwlabel{k};
    mask = zeros(size(im_tmp));
    tmp_lys = stack_lys{k};
    tmp_cal = stack_cal{k};
    tmp_bw = stack_bw{k};
    stats = regionprops(im_tmp,'all');
if k<10
    number=['00',num2str(k)];
elseif k<100
    number=['0',num2str(k)];
else
    number=num2str(k);
end
    

%% In each frame, only select largest object for quantification
obj_area = [];
for n = 1:size(stats,1)  
   obj_area = [obj_area,stats(n).Area]; 
end
obj_num = find(obj_area==max(obj_area));
all_obj_area = [all_obj_area,max(obj_area)];
all_obj_lys_intensity = [all_obj_lys_intensity,mean(tmp_lys(stats(obj_num).PixelIdxList))];

I = stack_bw{k};
% Measure shape parameters
boundaryPoint = 4;         %number of boundary points curvature is found over 
curvatureThresh = 0.8;     %the maximum allowed value of the curvature measure
bp_tangent = 3;            % number of boundary points the tangent angle is found over 
interpdmin = 0.3;           % the minimum number of pixels seperating boundary points after interpolation
loopclose = 1;              % 0 - if open boundaries | 1 - if closed boundaries

%% Find the curvature
[shape_details, Icurv] = curvature(I, obj_num, boundaryPoint, curvatureThresh, bp_tangent, ...
                                 interpdmin, loopclose);
        
%% Plot curvature
X = shape_details.XY(:,1);
Y = shape_details.XY(:,2);
Z = zeros(size(X));
C = shape_details.curvature'*1;
  

%% Calculate local intensity ratio value of a 3x3 region centered on boundary sample points
loc_ratio = [];
loc_cal = [];
mask = zeros(size(tmp_lys));
X_new=[];
Y_new=[];
C_new=[];
for m = 1:length(X)
    X0 = round(X(m)); Y0 = round(Y(m));
    loc_int_lys = [];
    loc_int_cal = [];
    if X0<=1 | Y0<=1 | X0==im_width | Y0==im_height
        continue;
    else
    for m1 = X0-thickness:X0+thickness
        if m1<1 | m1>im_width
             continue;
        else
        for m2 = Y0-thickness:Y0+thickness
            if m2<1 | m2>im_height
                continue;
            elseif tmp_bw(m2,m1)
            loc_int_lys = [loc_int_lys,tmp_lys(m2,m1)];
            loc_int_cal = [loc_int_cal,tmp_cal(m2,m1)];
            mask(m2,m1) = 100;
            end
        end
        end
    end
    loc_ratio = [loc_ratio,mean(loc_int_cal)/mean(loc_int_lys)];
    loc_cal = [loc_cal,mean(loc_int_cal)];
    X_new = [X_new;X(m)];
    Y_new = [Y_new;Y(m)];
    C_new = [C_new;C(m)]; %new curvature excluding image boundary points
    end
end

C2 = loc_ratio';

Z=zeros(size(X_new));
fig1 = figure('visible','off'); % changing to 'on' will somehow block the program 
imshow(Icurv);
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
hold on;
surf([X_new(:) X_new(:)], [Y_new(:) Y_new(:)], [Z(:) Z(:)], [C_new C_new], ...  % Reshape and replicate data
 'FaceColor', 'none', ...    % Don't bother filling faces with color
 'EdgeColor', 'interp', ...  % Use interpolated color for edges
 'LineWidth', color_width);            % Make a thicker line
hold off
cmap = jet;
colormap(cmap);
cb = colorbar;  % Add a colorbar
cb.Label.String = 'Curvature';
set(colorbar,'visible','off')

caxis([-0.4,0.4]);
filename1 = [dir_write1,'curv_',number,'.png'];
exportgraphics(fig1,filename1,'Resolution',export_resolution);


pause(0.1);
fig2=figure('visible','off'); 
imshow(Icurv);
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
hold on
surf([X_new(:) X_new(:)], [Y_new(:) Y_new(:)], [Z Z], [C2 C2], ...  % Reshape and replicate data
 'FaceColor', 'none', ...    % Don't bother filling faces with color
 'EdgeColor', 'interp', ...  % Use interpolated color for edges
 'LineWidth', color_width);            % Make a thicker line
hold off

cmap = jet;
colormap(cmap);
cb = colorbar;  % Add a colorbar
cb.Label.String = 'cal/lys intensity ratio';
set(colorbar,'visible','off')

caxis([0 2]);
filename2 = [dir_write2,'IntRatio_',number,'.png'];
exportgraphics(fig2,filename2,'Resolution',export_resolution);




%%calculate correlation coefficient between local intensity and curvature
%%at each time point

R1 = corrcoef(C_new(find(C_new>=0)),C2(find(C_new>=0))); 
% only quantify region with curvature>=0

r1 = R1(1,2);  %correlation coefficient between local curvature and intensity ratio

R_all=[R_all,r1];
end
 disp([im_name_base,' average correlation coefficient']);disp(mean(R_all));
 disp([im_name_base,' average main object area']);disp(mean(all_obj_area));
 disp([im_name_base,' average main object mean intensity']);disp(mean(all_obj_lys_intensity));
%% End parameters
%--------------------------------------------------------------------------
Runtime = toc(Start);

