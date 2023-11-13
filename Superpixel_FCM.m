%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Joined spatial and spectral segmentation of %%%%%%%
%%%% hyperspectral datasets on historical art objects %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%All the Matlab toolbox involved and need to be installed:
    %Image Processing Toolbox
    %Statistics and Machine Learning Toolbox
    %Image Processing Toolbox Hyperspectral Imaging Library 
    %from Add-on explorer
%Download from https://github.com/JakobSig/HSI2RGB : 
    %HSI2RGB.m and D_illuminants.mat 
    %put in the same working directory


%%%%%%%%%%%%%%%Pre-processing%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Two options : 1) read HSI cube
% DataHSI = multibandread('C:\Users\?\Documents\unsupervised_clustering\mock-up-Signac-training.bip', [599,456,172], 'float', 0, 'bip', 'ieee-le');
% DataLAB = multibandread('C:\Users\?\Documents\unsupervised_clustering\mock-up-Signac-training-CIE.bip', [599,456,3], 'float', 0, 'bip', 'ieee-le');
% enter in the brackets : [lines, samples, bands]

%Transform the LAB data into RGB data
% DataRGB = lab2rgb((DataLAB), WhitePoint="d65");%can choose the illuminant

% 2) or directly load the matrix
load('DataHSI.mat');
load('DataLAB.mat');
load('DataRGB.mat');

%read basic information of the datacube
[rows, cols, bands] = size(DataHSI);

%%%%%%%%%%%%%%superpixels%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define Number of Superpixel
NumSuperpixel = 580;
%run on 2D image-load the LAB or greyscale image
[LLAB,NLAB] = superpixels(DataLAB,NumSuperpixel,'IsInputLab',true,'Compactness',2);
% %LLAB: the superpixel label %NLAB: total number of superpixel generated 
% it is possible to modify the compactness

BWLAB = boundarymask(LLAB); %get boundary
Border = imoverlay(DataRGB,BWLAB,'cyan'); %image with boundary
figure;imshow(Border);
imwrite(Border,"superpixelborder580.png"); %save original resolution image

% Two ways to extract the centroid of each superpixel :

% % 1) spectral centroid: averaging within SP
% Centroid = zeros(NLAB,143);
% for j = 1:NLAB
%     [row, col] = find(LLAB == j);%extract spectrum indexes of the jth superpixel
%     C_SP = zeros(length(col),143);
%     for i = 1:length(col)
%         C_SP(i,:) = Signac(row(i),col(i),:);%extract the spectra
%     end
%     Centroid(j,:) = mean(C_SP); %average
% end

%spatial-centroids, take the spectrum in the center of the superpixel 
stats = regionprops(LLAB, 'Centroid');%get spatial center
index = cat(1, stats.Centroid);%extract the coordinates
idx_round = round(index);%round to integar

C_spatial = zeros(NLAB,bands);
for i = 1:NLAB
    C_spatial(i,:) = DataHSI(idx_round(i,2),idx_round(i,1),:);%extract spec at center position
end


%%%%%%%%%%%%%%%clustering%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%fuzzy clustering%%%%%
% Two ways to fixe the number of cluster :
% 1) Fixed by us :
NumClusters = 20;
% 2) Estimate by this function
% [centersN,sigmaN] = subclust(C_spatial,0.9);
% NumClusters = size(centersN,1);
options = [2 100 0.01 0]; %exponent,maxIteration,minImprovement,printSteps
[centers,U] = fcm(C_spatial,NumClusters,options);%cluster centers and possibility matrix
[maxU, clusters] = max(U);%assign clusters based on maximum

%label back to entire image
seg20 = LLAB;
for i = 1:NLAB
    idx = LLAB == i;
    seg20(idx) = clusters(i);%substitute superpixel belongingness to cluster membership
end

%load the document wavelength.mat :
load('wavelength.mat')

% obtain colormapValues 
[colormapValues,XYZ]=HSItoRGB(wavelength,centers,NumClusters,1,65,0.01);
colormapValues = reshape(colormapValues,NumClusters,3);
% plot the legend
% Function HSItoRGB is asked
figure; 
hold on;
for i = 1:NumClusters %number of labels
    rectangle('Position', [0, i, 1, 1], 'FaceColor', colormapValues(i, :));
    text(1.5, i + 0.5, sprintf('Cluster %d', i));
end
hold off;
ax = gca;
ax.Color = 'none'; % Set background color of the legend to transparent
ax.XLim = [0, 18];
ax.YLim = [0, NumClusters + 1]; %change with number of labels
ax.Visible = 'off';

%plot cluster centers
figure;
hold on;
for i = 1:NumClusters
    plot(wavelength, centers(i, :), 'Color', colormapValues(i, :),'LineWidth',1.5);
end
hold off;

xlabel('Wavelength');
ylabel('Reflectance');
title('Cluster centers');
legend('Cluster 1', 'Cluster 2', 'Cluster 3', ...
    'Cluster 4','Cluster 5','Cluster 6','Cluster 7', ...
    'Cluster 8','Cluster 9','Cluster 10', ...
    'Cluster 11','Cluster 12','Cluster 13','Cluster 14','Cluster 15','Cluster 16','Cluster 17','Cluster 18','Clsuter 19','Cluster 20');
% save spectra

%generate the maps
%all clusters
labelRGB = label2rgb(seg20, colormapValues);
figure;imshow(labelRGB)
imwrite(labelRGB,"labelmapFCM.png");

%specific cluster
% white = ones(rows,cols,3)*255;
mapi = labeloverlay(DataRGB,seg20,"Colormap",colormapValues,"IncludedLabels",[4,5],"Transparency",0);
figure;imshow(mapi)
imwrite(mapi,"clusters4_5_FCM.png");
%can replace DataRGB by white, to have a white background

%automated process for each single cluster
% for i = 1:NumClusters
%     mapi = labeloverlay(DataRGB,seg20,"Colormap",colormapValues,"IncludedLabels",[i],"Transparency",0);
%     fileName = sprintf('cluster%d_FCM.png', i);
%     imwrite(mapi, fileName);
% end